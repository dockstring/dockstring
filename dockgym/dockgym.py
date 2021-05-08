import os
import platform
import subprocess
import sys
import tempfile
from pathlib import Path

import pkg_resources
from rdkit import rdBase
from rdkit.Chem import AllChem as Chem

from .utils import DockingError


def load_target(name):
    return Target(name)


class Target:
    def __init__(self, name):
        self.name = name

        self._bin_dir = Path(pkg_resources.resource_filename(__package__, 'bin')).resolve()
        self._receptors_dir = Path(pkg_resources.resource_filename(__package__, 'receptors')).resolve()

        self._vina = self._bin_dir / self._get_vina_filename()

        # Create temporary directory where the PDB, PDBQT and conf files for the target will be saved
        self._tmp_dir_handle = tempfile.TemporaryDirectory()
        self._tmp_dir = Path(self._tmp_dir_handle.name).resolve()

        # Set PDB, PDBQT, and conf files
        self._pdb = self._receptors_dir / (self.name + '_receptor.pdb')
        self._pdbqt = self._receptors_dir / (self.name + '_receptor.pdbqt')
        self._conf = self._receptors_dir / (self.name + '_conf.txt')

        # Ensure files exist
        if not all(p.exists() for p in [self._tmp_dir, self._pdb, self._pdbqt, self._conf]):
            raise DockingError(f"'{self.name}' is not a target we support")

    @staticmethod
    def _get_vina_filename() -> str:
        system_name = platform.system()
        if system_name == 'Linux':
            return 'vina_linux'
        else:
            raise DockingError(f"System '{system_name}' not yet supported")

    @staticmethod
    def _smiles_or_inchi_to_mol(smiles_or_inchi, verbose=False):
        if not verbose:
            rdBase.DisableLog('rdApp.error')

        mol = Chem.MolFromSmiles(smiles_or_inchi)
        if mol is None:
            mol = Chem.MolFromInchi(smiles_or_inchi)
            if mol is None:
                raise DockingError('Could not parse SMILES or InChI string')

        if not verbose:
            rdBase.EnableLog('rdApp.error')

        return mol

    @staticmethod
    def _embed_mol(mol, seed: int, max_num_attempts: int = 10):
        """Will attempt to find 3D coordinates <max_num_attempts> times with different random seeds"""

        # Add hydrogen atoms in order to get a sensible 3D structure, and remove them later
        mol = Chem.AddHs(mol)
        Chem.EmbedMolecule(mol, randomSeed=seed, maxAttempts=max_num_attempts)
        mol = Chem.RemoveHs(mol)  # TODO: Why?

        # If not a single conformation is obtained in all the attempts, raise an error
        if len(mol.GetConformers()) == 0:
            raise DockingError('Generation of ligand conformation failed')

        return mol

    @staticmethod
    def _write_embedded_mol_to_pdb(mol, ligand_pdb):
        if len(mol.GetConformers()) < 1:
            raise DockingError('For conversion to PDB a conformer is required')
        Chem.MolToPDBFile(mol, filename=str(ligand_pdb))

    @staticmethod
    def _convert_pdb_to_pdbqt(ligand_pdb, ligand_pdbqt, verbose=False):
        # yapf: disable
        cmd_list = [
            'obabel',
            '-ipdb', ligand_pdb,
            '-opdbqt',
            '-O', ligand_pdbqt,
            '--partialcharge', 'gasteiger'
        ]
        # yapf: enable
        cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = cmd_return.stdout.decode('utf-8')

        # If verbose, print output to string
        if verbose:
            print(output)

        # If failure, raise DockingError
        if cmd_return.returncode != 0:
            raise DockingError('Conversion from PDB to PDBQT failed')

    def _dock_pdbqt(self, ligand_pdbqt, vina_logfile, vina_outfile, seed, num_cpu=1, verbose=False):
        # yapf: disable
        cmd_list = [
            str(self._vina),
            '--receptor', self._pdbqt,
            '--config', self._conf,
            '--ligand', ligand_pdbqt,
            '--log', vina_logfile,
            '--out', vina_outfile,
            '--seed', str(seed),
        ]
        # yapf: enable
        if num_cpu is not None:
            cmd_list += ['--cpu', str(num_cpu)]

        cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = cmd_return.stdout.decode('utf-8')

        # If verbose, print output to string
        if verbose:
            print(output)

        # If failure, raise DockingError
        if cmd_return.returncode != 0:
            raise DockingError('Docking with Vina failed')

    @staticmethod
    def _get_top_score_from_vina_logfile(vina_logfile):
        try:
            with open(vina_logfile, mode='r') as f:
                counter_to_score = None
                for each_line in f:
                    # Try to find the table header. Once found, count three lines to get the score
                    if 'mode |   affinity | dist from best mode' in each_line:
                        counter_to_score = 0
                    elif counter_to_score == 2:
                        line_with_score = each_line
                        break
                    elif counter_to_score is not None:
                        counter_to_score += 1
            return float(line_with_score.split()[1])

        except Exception:
            raise DockingError('Failed to find suitable poses')

    def dock(self, mol, num_cpu=1, seed=974528263, verbose=False):
        """
        Given a molecule, this method will return a docking score against the current target.
        - mol: either a SMILES string, an InChIKey or a RDKit molecule object
        - num_cpu: number of cpus that AutoDock Vina should use for the docking. By default,
          it will try to find all the cpus on the system, and failing that, it will use 1.
        - seed: integer random seed for reproducibility

        The process is the following:
        1. Obtain RDKit molecule object
        2. Embed molecule to 3D conformation
        3. Prepare ligand
        4. Dock
        5. Extract all the info from the docking output
        """

        # Docking with Vina is performed in a temporary directory
        ligand_pdb = self._tmp_dir / 'ligand.pdb'
        ligand_pdbqt = self._tmp_dir / 'ligand.pdbqt'
        vina_logfile = self._tmp_dir / 'vina.log'
        vina_outfile = self._tmp_dir / 'vina.out'

        try:
            # Prepare ligand
            if not isinstance(mol, Chem.Mol):
                mol = self._smiles_or_inchi_to_mol(mol, verbose=verbose)
                mol = self._embed_mol(mol, seed=seed)

            # Prepare ligand files
            self._write_embedded_mol_to_pdb(mol, ligand_pdb)
            self._convert_pdb_to_pdbqt(ligand_pdb, ligand_pdbqt, verbose=verbose)

            # Dock
            self._dock_pdbqt(ligand_pdbqt, vina_logfile, vina_outfile, seed=seed, num_cpu=num_cpu, verbose=verbose)

            # Process docking output
            score = self._get_top_score_from_vina_logfile(vina_logfile)

            # TODO Get the pose from vina_outfile

            return score, None

        except DockingError as error:
            mol_id = Chem.MolToSmiles(mol) if isinstance(mol, Chem.Mol) else mol
            raise DockingError(f"An error occurred for ligand '{mol_id}': {error}")

        # TODO Include Mac and Windows binaries in the repository
        # TODO Put all the calculated scores (and maybe the poses too?) under "data". What should be the format?
        # - Plain text for the scores and smiles?
        # - What format for the poses?

    def info(self):
        """
        Print some info about the target.
        """
        pass

    def view(self, search_box=True):
        """
        Start pymol and view the receptor and the search box.
        """
        with open(self._conf, 'r') as f:
            # Extract the search box information
            for line in f:
                if 'center_x' in line:
                    center_x = float(line.split()[2])
                elif 'center_y' in line:
                    center_y = float(line.split()[2])
                elif 'center_z' in line:
                    center_z = float(line.split()[2])
                elif 'size_x' in line:
                    size_x = float(line.split()[2])
                elif 'size_y' in line:
                    size_y = float(line.split()[2])
                elif 'size_z' in line:
                    size_z = float(line.split()[2])

            receptor = str(self._pdb.resolve())
            # TODO Remove ligand from the main view() method, make it a part of the docked ligands' method
            ligand = f'/home/mgarort/repos/dockgym/playground/crystal_ligands/{self.name}/crystal_ligand.mol2'
            pymol_view_search_box = Path(
                sys.modules[self.__class__.__module__].__file__).parent.parent / 'utils' / 'view_search_box.py'
            command = f'pymol -R {pymol_view_search_box} {receptor} {ligand}'
            if search_box:
                command += f" -d 'view_search_box center_x={center_x}, center_y={center_y}, center_z={center_z}, \
                              size_x={size_x}, size_y={size_y}, size_z={size_z}'"

            command += ' > /dev/null'
            # TODO Change to subprocess
            os.system(command)
