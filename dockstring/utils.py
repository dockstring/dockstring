import logging
import os
import platform
import re
import subprocess
from pathlib import Path
from typing import List, Union, Dict

import pkg_resources
from rdkit import rdBase
from rdkit.Chem import AllChem as Chem

PathType = Union[str, os.PathLike]


class DockingError(Exception):
    """Raised when Target.dock fails at any step"""
    pass


def get_vina_filename() -> str:
    system_name = platform.system()
    if system_name == 'Linux':
        return 'vina_linux'
    else:
        raise DockingError(f"System '{system_name}' not yet supported")


def get_resources_dir() -> Path:
    path = Path(pkg_resources.resource_filename(__package__, 'resources'))
    if not path.is_dir():
        raise DockingError("'resources' directory not found")
    return path


def get_targets_dir() -> Path:
    path = get_resources_dir() / 'targets'
    if not path.is_dir():
        raise DockingError("'targets' directory not found")
    return path


def get_bin_dir() -> Path:
    path = get_resources_dir() / 'bin'
    if not path.is_dir():
        raise DockingError("'bin' directory not found")
    return path


def get_vina_path() -> Path:
    path = get_bin_dir() / get_vina_filename()
    if not path.is_file():
        raise DockingError('AutoDock Vina executable not found')
    return path


def smiles_or_inchi_to_mol(smiles_or_inchi, verbose=False):
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


def embed_mol(mol, seed: int, max_num_attempts: int = 10):
    """Will attempt to find 3D coordinates <max_num_attempts> times with different random seeds"""
    # Add hydrogen atoms in order to get a sensible 3D structure
    mol = Chem.AddHs(mol)
    Chem.EmbedMolecule(mol, randomSeed=seed, maxAttempts=max_num_attempts)
    # If not a single conformation is obtained in all the attempts, raise an error
    if mol.GetNumConformers() == 0:
        raise DockingError('Generation of ligand conformation failed')
    return mol


def refine_mol_with_ff(mol, max_iters=1000):
    """
    Will attempt to refine the embedded coordinates.
    """
    try:
        opt_failed = Chem.MMFFOptimizeMolecule(mol, mmffVariant='MMFF94', maxIters=max_iters)
    except Chem.rdchem.KekulizeException as exception:
        raise DockingError('Structure refinement of ligand failed because the ligand could not be kekulized.\n'
                           'Message by RDKit:\n'
                           f'{exception}')
    if opt_failed != 0:
        raise DockingError('Structure refinement of ligand failed')


def write_embedded_mol_to_pdb(mol, ligand_pdb):
    if mol.GetNumConformers() < 1:
        raise DockingError('For conversion to PDB a conformer is required')
    Chem.MolToPDBFile(mol, filename=str(ligand_pdb))


def convert_pdbqt_to_pdb(pdbqt_file: PathType, pdb_file: PathType, verbose=False) -> None:
    # yapf: disable
    cmd_args = [
        'obabel',
        '-ipdbqt', pdbqt_file,
        '-opdb',
        '-O', pdb_file,
    ]
    # yapf: enable
    cmd_return = subprocess.run(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout = cmd_return.stdout.decode('utf-8')

    if verbose:
        logging.info(stdout)

    if cmd_return.returncode != 0:
        raise DockingError('Conversion from PDBQT to PDB failed')


def protonate_pdb(pdb_file: PathType, verbose=False):
    # Remove all hydrogen atoms
    # yapf: disable
    cmd_list = [
        'obabel',
        pdb_file,
        '-O', pdb_file,
        '-d',  # delete hydrogen atoms
    ]
    # yapf: enable
    cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = cmd_return.stdout.decode('utf-8')

    if verbose:
        logging.info(output)

    if cmd_return.returncode != 0:
        raise DockingError('Protonation of ligand failed.')

    # Add hydrogen atoms for pH 7
    # yapf: disable
    cmd_list = [
        'obabel',
        pdb_file,
        '-O', pdb_file,
        '-p', '7.0',  # add hydrogen atoms
    ]
    # yapf: enable
    cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = cmd_return.stdout.decode('utf-8')

    if verbose:
        logging.info(output)

    if cmd_return.returncode != 0:
        raise DockingError('Protonation of ligand failed.')


def convert_pdb_to_pdbqt(pdb_file: PathType, pdbqt_file: PathType, verbose=False):
    # yapf: disable
    cmd_list = [
        'obabel',
        '-ipdb', pdb_file,
        '-opdbqt',
        '-O', pdbqt_file,
        '--partialcharge', 'gasteiger'
    ]
    # yapf: enable
    cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = cmd_return.stdout.decode('utf-8')

    if verbose:
        logging.info(output)

    if cmd_return.returncode != 0:
        raise DockingError('Conversion from PDB to PDBQT failed')


def read_mol_from_pdb(pdb_file: PathType) -> Chem.Mol:
    return Chem.MolFromPDBFile(str(pdb_file))


def write_mol_to_pdb(mol: Chem.Mol, pdb_file: PathType):
    if mol.GetNumConformers() < 1:
        raise DockingError('For conversion to PDB a conformer is required')
    Chem.MolToPDBFile(mol, filename=str(pdb_file))


real_number_pattern = r'[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?'
score_re = re.compile(rf'REMARK VINA RESULT:\s*(?P<score>{real_number_pattern})')


def parse_scores_from_pdb(pdb_file: PathType) -> List[float]:
    with open(pdb_file, mode='r') as f:
        content = f.read()
    return [float(match.group('score')) for match in score_re.finditer(content)]


conf_re = re.compile(rf'^(?P<key>\w+)\s*=\s*(?P<value>{real_number_pattern})\s*\n$')


def parse_search_box_conf(conf_file: PathType) -> Dict[str, float]:
    d = {}
    with open(conf_file, mode='r') as f:
        for line in f.readlines():
            match = conf_re.match(line)
            if match:
                d[match.group('key')] = float(match.group('value'))

        assert len(d) == 6
        return d
