import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pkg_resources
from rdkit import rdBase
from rdkit.Chem import AllChem as Chem

from .utils import DockingError, get_num_conf


def load_target(name):
    return Target(name)


class Target:
    def __init__(self, name, random_seed=974528263):
        self.name = name
        self.random_seed = random_seed
        # TODO Determine whether the OS is Linux, Mac or Windows
        # # Define paths to dependencies
        # system = platform.system()
        # if system == 'Linux':
        #     # vina_binary = lib / 'vina_linux'
        #     # pythonsh_binary = lib / 'pythonsh'
        #     # TODO Check if pythonsh is platform-dependent, or it is the same for Linux, Mac, Windows, etc
        #     #      I can do this with `diff`
        #     pass
        # elif system == 'Windows':
        #     pass
        # elif system == 'Darwin':
        #     pass

        self._bin_dir = Path(pkg_resources.resource_filename(__package__, 'bin')).resolve()
        self._receptors_dir = Path(pkg_resources.resource_filename(__package__, 'receptors')).resolve()

        self._vina = self._bin_dir / 'vina_linux'

        # Create temporary directory where the PDB, PDBQT and conf files for the target will be saved
        self._tmp_dir_handle = tempfile.TemporaryDirectory()
        self._tmp_dir = Path(self._tmp_dir_handle.name).resolve()

        # Copy PDB, PDBQT and conf files to the temporary directory
        # TODO Create PDBQT for receptors
        self._pdb = self._receptors_dir / (self.name + '_receptor.pdb')
        self._pdbqt = self._receptors_dir / (self.name + '_receptor.pdbqt')
        self._conf = self._receptors_dir / (self.name + '_conf.txt')

    def __del__(self):
        self._tmp_dir_handle.cleanup()

    # Submethods for docking a ligand
    def _smiles_or_inchi_2_mol(self, smiles_or_inchi):
        rdBase.DisableLog('rdApp.error')
        mol = Chem.MolFromSmiles(smiles_or_inchi)
        if mol is None:
            mol = Chem.MolFromInchi(smiles_or_inchi)
            if mol is None:
                raise DockingError(
                    f'Docking of molecule {self._mol_id} failed because '
                    'it is not a proper SMILES, InChi or RDKit molecule object.'
                )
        rdBase.EnableLog('rdApp.error')
        return mol

    def _mol_2_embedding(self, mol):
        # Will attempt to find 3D coordinates 10 times with different random seeds
        n_attempts = 10
        # Set first random seed for reproducibility
        random_seed = self._dock_random_seed
        for i in range(n_attempts):
            # Simple approach to get subsequent random seeds from first random seed
            random_seed = (random_seed -
                           (-1)**(i % 2) * i * random_seed) % 766523564
            # Always add hydrogens in order to get a sensible 3D structure, and remove them later
            mol = Chem.AddHs(mol)
            Chem.EmbedMolecule(mol, randomSeed=random_seed)
            mol = Chem.RemoveHs(mol)
            # If at least one conformation has been obtained, don't try more and break out of the loop.
            # Otherwise, keep trying to hopefully generate a ligand conformation
            if get_num_conf(mol) > 0:
                break
        # If not a single conformation is obtained in all the attempts, raise an error and return None.
        # Otherwise, return the molecule with the conformation
        if get_num_conf(mol) == 0:
            raise DockingError(
                f'Docking of molecule  {self._mol_id}  failed during the '
                'ligand conformation generation with RDKit.')
        else:
            return mol

    def _embedding_2_pdb(self, mol, ligand_pdb):
        Chem.MolToPDBFile(mol, str(ligand_pdb))

    def _pdb_to_pdbqt(self, ligand_pdb, ligand_pdbqt):
        # yapf: disable
        cmd_list = [
            'obabel',
            '-ipdb', ligand_pdb,
            '-opdbqt',
            '-O', ligand_pdbqt,
            '--partialcharge', 'gasteiger'
        ]
        # yapf: enable
        cmd_return = subprocess.run(cmd_list,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT)
        output = cmd_return.stdout.decode('utf-8')

        # If there exists a logfile handle, save output there. If not, do nothing
        try:
            self._dock_logfile_handle.write(output)
        except AttributeError as exception:
            pass

        # If verbose, print output to string
        if self._dock_verbose:
            print(output)

        # If failure, raise DockingError
        if cmd_return.returncode != 0:
            raise DockingError(f'Docking of molecule {self._mol_id} failed '
                               f'during the conversion of PDB to PDBQT with OpenBabel.')

    def _dock_pdbqt(self,
                    ligand_pdbqt,
                    vina_logfile,
                    vina_outfile,
                    num_cpu=None):
        # yapf: disable
        cmd_list = [
            str(self._vina),
            '--receptor', self._pdbqt,
            '--config', self._conf,
            '--ligand', ligand_pdbqt,
            '--log', vina_logfile,
            '--out', vina_outfile,
            '--seed', str(self._dock_random_seed),
        ]
        # yapf: enable
        if num_cpu is not None:
            cmd_list += ['--cpu', str(num_cpu)]

        cmd_return = subprocess.run(cmd_list,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT)
        output = cmd_return.stdout.decode('utf-8')

        # If there exists a logfile handle, save output there. If not, do nothing
        try:
            self._dock_logfile_handle.write(output)
        except AttributeError as exception:
            pass

        # If verbose, print output to string
        if self._dock_verbose:
            print(output)

        # If failure, raise DockingError
        if cmd_return.returncode != 0:
            raise DockingError(
                f'Docking of molecule {self._mol_id} failed during the docking with Vina.'
            )

    def _get_top_score_from_vina_logfile(self, vina_logfile):
        try:
            with open(vina_logfile, 'r') as f:
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
            top_score = float(line_with_score.split()[1])
            return top_score
        except Exception as exception:
            raise DockingError(f'Docking of molecule {self._mol_id}  failed '
                               'because no suitable poses were found.')

    def dock(self, mol, num_cpu=None, seed=None, logfile=None, verbose=False):
        """
        Given a molecule, this method will return a docking score against the current target.
        - mol: either a SMILES string, an inchikey or a RDKit molecule object
        - num_cpu: number of cpus that Autodock Vina should use for the docking. By default,
          it will try to find all the cpus on the system, and failing that, it will use 1.
        - seed: integer random seed for reproducibility

        The process is the following:
        1. Obtain RDKit molecule object
        2. Embed molecule to 3D conformation
        3. Prepare ligand
        4. Dock
        5. Extract all the info from the docking output
        """
        # Define molecule identifier for error messages
        if isinstance(mol, Chem.Mol):
            self._mol_id = Chem.MolToSmiles(mol)
        else:
            self._mol_id = mol
        if seed is None:
            self._dock_random_seed = self.random_seed
        else:
            self._dock_random_seed = seed
        if logfile is not None:
            self._dock_logfile_handle = open(logfile, 'w')
        self._dock_verbose = verbose
        # Docking with Vina is performed in a temporary directory
        with tempfile.TemporaryDirectory() as dock_tmp_dir:
            # Define necessary filenames
            dock_tmp_dir = Path(dock_tmp_dir)
            ligand_pdb = (dock_tmp_dir / 'ligand.pdb').resolve()
            ligand_pdbqt = (dock_tmp_dir / 'ligand.pdbqt').resolve()
            vina_logfile = (dock_tmp_dir / 'vina.log').resolve()
            vina_outfile = (dock_tmp_dir / 'vina.out').resolve()
            try:
                # Prepare ligand
                # TODO Handle RDKit output too with verbose/logfile
                if not isinstance(mol, Chem.Mol):
                    mol = self._smiles_or_inchi_2_mol(mol)
                mol = self._mol_2_embedding(mol)
                self._embedding_2_pdb(mol, ligand_pdb)
                self._pdb_to_pdbqt(ligand_pdb, ligand_pdbqt)
                # Dock
                self._dock_pdbqt(ligand_pdbqt,
                                 vina_logfile,
                                 vina_outfile,
                                 num_cpu=num_cpu)
                # Process docking output
                score = self._get_top_score_from_vina_logfile(vina_logfile)
                # TODO Get the pose from vina_outfile

                # Clean up and return
                del self._dock_random_seed
                if hasattr(self, '_dock_logfile_handle'):
                    self._dock_logfile_handle.close()
                    del self._dock_logfile_handle
                del self._dock_verbose
                return score, None
            except DockingError as error:
                print(f'DockingError: ' + str(error))
                del self._dock_random_seed
                if hasattr(self, '_dock_logfile_handle'):
                    self._dock_logfile_handle.close()
                    del self._dock_logfile_handle
                del self._dock_verbose
                return None, None

            # Using tempfile makes the temporary file creation portable for all systems
            print(dock_tmp_dir)

        # TODO Include Mac and Windows binaries in the repository
        # TODO Check whether the current system is Linux, Mac or Windows with platform.system()
        #      https://stackoverflow.com/questions/22321397/python-os-name-return-nt-on-windows-7
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
                sys.modules[self.__class__.__module__].__file__
            ).parent.parent / 'utils' / 'view_search_box.py'
            command = f'pymol -R {pymol_view_search_box} {receptor} {ligand}'
            if search_box:
                command += f" -d 'view_search_box center_x={center_x}, center_y={center_y}, center_z={center_z}, \
                              size_x={size_x}, size_y={size_y}, size_z={size_z}'"

            command += ' > /dev/null'
            # TODO Change to subprocess
            os.system(command)
