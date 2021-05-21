import logging
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, List

from rdkit.Chem import AllChem as Chem

from .utils import (DockingError, smiles_or_inchi_to_mol, embed_mol, refine_mol_with_ff, write_embedded_mol_to_pdb,
                    protonate_pdb, convert_pdbqt_to_pdb, convert_pdb_to_pdbqt, read_mol_from_pdb, parse_scores_from_pdb,
                    parse_search_box_conf, PathType, get_targets_dir, get_vina_path, get_resources_dir)

logging.basicConfig(format='%(message)s')


def load_target(name: str, *args, **kwargs):
    return Target(name, *args, **kwargs)


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
    def __init__(self, name, working_dir: Optional[PathType] = None):
        self.name = name

        # Directory where the ligand and output files will be saved
        self._custom_working_dir = working_dir
        self._tmp_dir_handle: Optional[tempfile.TemporaryDirectory] = None

        # Set PDB, PDBQT, and conf files
        targets_dir = get_targets_dir()
        self._pdb = targets_dir / (self.name + '_target.pdb')
        self._pdbqt = targets_dir / (self.name + '_target.pdbqt')
        self._conf = targets_dir / (self.name + '_conf.txt')

        # Ensure files exist
        if not all(p.exists() for p in [self._pdb, self._pdbqt, self._conf]):
            raise DockingError(f"'{self.name}' is not a target we support")

    def __repr__(self):
        return f"dockstring.Target(name='{self.name}', dir='{self._tmp_dir}')"

    @property
    def _tmp_dir(self) -> Path:
        if self._custom_working_dir:
            return Path(self._custom_working_dir).resolve()

        # If no custom working dir is set and the tmp working dir handle is not initialized, initialize it
        if not self._tmp_dir_handle:
            self._tmp_dir_handle = tempfile.TemporaryDirectory()

        return Path(self._tmp_dir_handle.name).resolve()

    def _dock_pdbqt(self,
                    ligand_pdbqt,
                    vina_logfile,
                    vina_outfile,
                    seed,
                    num_cpus: Optional[int] = None,
                    verbose=False):
        # yapf: disable
        cmd_list = [
            get_vina_path(),
            '--receptor', self._pdbqt,
            '--config', self._conf,
            '--ligand', ligand_pdbqt,
            '--log', vina_logfile,
            '--out', vina_outfile,
            '--seed', str(seed),
        ]
        # yapf: enable
        if num_cpus is not None:
            cmd_list += ['--cpu', str(num_cpus)]

        cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = cmd_return.stdout.decode('utf-8')

        if verbose:
            logging.info(output)

        # If failure, raise DockingError
        if cmd_return.returncode != 0:
            raise DockingError('Docking with Vina failed')

    def dock(self, string: str, num_cpus: Optional[int] = None, seed=974528263, verbose=False):
        """
        Given a molecule, this method will return a docking score against the current target.
        - mol: either a SMILES or an InChI string
        - num_cpus: number of cpus that AutoDock Vina should use for the docking. By default,
          it will try to find all the cpus on the system, and failing that, it will use 1.
        - seed: integer random seed for reproducibility

        The process is the following:
        1. Obtain RDKit molecule object
        2. Embed molecule to 3D conformation
        3. Refine embedding with a force field.
        4. Prepare ligand
        6. Dock
        8. Parse docking results from the output
        """

        # Docking with Vina is performed in a temporary directory
        ligand_pdb = self._tmp_dir / 'ligand.pdb'
        ligand_pdbqt = self._tmp_dir / 'ligand.pdbqt'
        vina_logfile = self._tmp_dir / 'vina.log'
        vina_outfile = self._tmp_dir / 'vina.out'
        docked_ligand_pdb = self._tmp_dir / 'docked_ligand.pdb'

        try:
            # Prepare ligand
            mol = smiles_or_inchi_to_mol(string, verbose=verbose)
            embedded_mol = embed_mol(mol, seed=seed)
            refine_mol_with_ff(embedded_mol)

            # Prepare ligand files
            write_embedded_mol_to_pdb(embedded_mol, ligand_pdb)
            protonate_pdb(ligand_pdb, verbose=verbose)
            convert_pdb_to_pdbqt(ligand_pdb, ligand_pdbqt, verbose=verbose)

            # Dock
            self._dock_pdbqt(ligand_pdbqt, vina_logfile, vina_outfile, seed=seed, num_cpus=num_cpus, verbose=verbose)

            # Process docking output
            # If Vina does not find any appropriate poses, the output file will be empty
            if os.stat(vina_outfile).st_size == 0:
                raise DockingError('AutoDock Vina could not find any appropriate pose.')

            convert_pdbqt_to_pdb(pdbqt_file=vina_outfile, pdb_file=docked_ligand_pdb, verbose=verbose)
            ligands = read_mol_from_pdb(docked_ligand_pdb)
            scores = parse_scores_from_pdb(docked_ligand_pdb)

            if ligands is not None:
                assert len(scores) == ligands.GetNumConformers()

            return scores[0], {
                'ligands': ligands,
                'scores': scores,
            }

        except DockingError as error:
            logging.error(f"An error occurred for ligand '{string}': {error}")
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

    def view(self, mols: List[Chem.Mol] = None, search_box=True):
        """
        Start pymol and view the receptor and the search box.
        """
        commands = ['pymol', self._pdb]

        if search_box:
            pymol_script = get_resources_dir() / 'view_search_box.py'
            conf = parse_search_box_conf(self._conf)
            # yapf: disable
            commands += [
                pymol_script,
                '-d', 'view_search_box center_x={center_x}, center_y={center_y}, center_z={center_z}, '
                      'size_x={size_x}, size_y={size_y}, size_z={size_z}'.format(**conf)
            ]
            # yapf: enable

        if mols:
            tmp_dir_handle = tempfile.TemporaryDirectory()
            tmp_dir = Path(tmp_dir_handle.name).resolve()

            for index, mol in enumerate(mols):
                mol_pdb_file = tmp_dir / f'ligand_{index}.pdb'
                write_embedded_mol_to_pdb(mol, mol_pdb_file)
                commands += [str(mol_pdb_file)]

        return subprocess.run(commands)
