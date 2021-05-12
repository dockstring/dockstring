import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, List

import pkg_resources
from rdkit.Chem import AllChem as Chem

from dockgym.utils import (DockingError, get_vina_filename, smiles_or_inchi_to_mol, embed_mol,
                           write_embedded_mol_to_pdb, convert_pdbqt_to_pdb, convert_pdb_to_pdbqt, read_pdb_to_mol,
                           parse_scores_from_pdb, parse_search_box_conf)

import logging
logging.basicConfig(format='%(message)s')

def get_targets_dir() -> Path:
    return Path(pkg_resources.resource_filename(__package__, 'targets')).resolve()


def load_target(name):
    return Target(name)


def list_all_target_names() -> List[str]:
    targets_dir = get_targets_dir()
    file_names = [f for f in os.listdir(targets_dir) if os.path.isfile(os.path.join(targets_dir, f))]

    target_re = re.compile(r'^(?P<name>\w+)_target\.pdb$')
    names = []
    for file_name in file_names:
        match = target_re.match(file_name)
        if match:
            names.append(match.group('name'))

    return names


class Target:
    def __init__(self, name):
        self.name = name

        self._bin_dir = Path(pkg_resources.resource_filename(__package__, 'bin')).resolve()
        self._targets_dir = get_targets_dir()

        self._vina = self._bin_dir / get_vina_filename()

        # Create temporary directory where the PDB, PDBQT and conf files for the target will be saved
        self._tmp_dir_handle: Optional[tempfile.TemporaryDirectory] = None

        # Set PDB, PDBQT, and conf files
        self._pdb = self._targets_dir / (self.name + '_target.pdb')
        self._pdbqt = self._targets_dir / (self.name + '_target.pdbqt')
        self._conf = self._targets_dir / (self.name + '_conf.txt')

        # Ensure files exist
        if not all(p.exists() for p in [self._pdb, self._pdbqt, self._conf]):
            raise DockingError(f"'{self.name}' is not a target we support")

    @property
    def _tmp_dir(self) -> Path:
        if not self._tmp_dir_handle:
            self._tmp_dir_handle = tempfile.TemporaryDirectory()
        return Path(self._tmp_dir_handle.name).resolve()

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

        if verbose:
            logging.info(output)

        # If failure, raise DockingError
        if cmd_return.returncode != 0:
            raise DockingError('Docking with Vina failed')

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
        docked_ligand_pdb = self._tmp_dir / 'docked_ligand.pdb'

        try:
            # Prepare ligand
            if not isinstance(mol, Chem.Mol):
                mol = smiles_or_inchi_to_mol(mol, verbose=verbose)
                mol = embed_mol(mol, seed=seed)

            # Prepare ligand files
            write_embedded_mol_to_pdb(mol, ligand_pdb)
            convert_pdb_to_pdbqt(ligand_pdb, ligand_pdbqt, verbose=verbose)

            # Dock
            self._dock_pdbqt(ligand_pdbqt, vina_logfile, vina_outfile, seed=seed, num_cpu=num_cpu, verbose=verbose)

            # Process docking output
            # If Vina does not find any appropriate poses, the output file will be empty
            if os.stat(vina_outfile).st_size == 0:
                raise DockingError('AutoDock Vina could not find any appropriate pose.')

            convert_pdbqt_to_pdb(pdbqt_file=vina_outfile, pdb_file=docked_ligand_pdb, verbose=verbose)
            ligands = read_pdb_to_mol(docked_ligand_pdb)
            scores = parse_scores_from_pdb(docked_ligand_pdb)

            assert len(scores) == len(ligands.GetConformers())

            return scores[0], {
                'ligands': ligands,
                'scores': scores,
            }

        except DockingError as error:
            mol_id = Chem.MolToSmiles(mol) if isinstance(mol, Chem.Mol) else mol
            logging.error(f"DockingError: An error occurred for ligand '{mol_id}': {error}")
            return (None, None)

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
        Start pymol and view the target and the search box.
        """
        pymol_view_search_box_file = pkg_resources.resource_filename(__package__,
                                                                     os.path.join('utils', 'view_search_box.py'))
        commands = ['pymol', pymol_view_search_box_file, self._pdb]

        if search_box:
            conf = parse_search_box_conf(self._conf)
            commands += [
                '-d', 'view_search_box center_x={center_x}, center_y={center_y}, center_z={center_z}, '
                'size_x={size_x}, size_y={size_y}, size_z={size_z}'.format(**conf)
            ]

        return subprocess.run(commands)
