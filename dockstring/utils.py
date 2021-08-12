import copy
import logging
import os
import platform
import re
import subprocess
import sys
from pathlib import Path
from typing import List, Union, Dict, Optional

import pkg_resources
from rdkit import rdBase
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.Descriptors import NumRadicalElectrons
from rdkit.Chem.MolStandardize.rdMolStandardize import Uncharger

from .errors import (DockstringError, CanonicalizationError, ParsingError, SanityError, EmbeddingError,
                     StructureOptimizationError, FormatConversionError, ProtonationError, DockingError,
                     PoseProcessingError, OutputError)

PathType = Union[str, os.PathLike]


def setup_logger(level: Union[int, str] = logging.INFO, path: Optional[str] = None):
    logger = logging.getLogger("dockstring")
    logger.setLevel(level)

    formatter = logging.Formatter('%(asctime)s.%(msecs)03d %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if path is not None:
        fh = logging.FileHandler(path)
        fh.setFormatter(formatter)
        logger.addHandler(fh)


def get_vina_filename() -> str:
    system_name = platform.system()
    if system_name == 'Linux':
        return 'vina_linux'
    else:
        raise DockstringError(f"System '{system_name}' not yet supported")


def get_resources_dir() -> Path:
    path = Path(pkg_resources.resource_filename(__package__, 'resources'))
    if not path.is_dir():
        raise DockstringError("'resources' directory not found")
    return path


def get_targets_dir() -> Path:
    path = get_resources_dir() / 'targets'
    if not path.is_dir():
        raise DockstringError("'targets' directory not found")
    return path


def get_bin_dir() -> Path:
    path = get_resources_dir() / 'bin'
    if not path.is_dir():
        raise DockstringError("'bin' directory not found")
    return path


def get_vina_path() -> Path:
    path = get_bin_dir() / get_vina_filename()
    if not path.is_file():
        raise DockstringError('AutoDock Vina executable not found')
    return path


def canonicalize_smiles(smiles: str) -> str:
    try:
        return Chem.CanonSmiles(smiles, useChiral=True)
    except Exception as e:
        raise CanonicalizationError(f'Cannot canonicalize SMILES: {e}')


def smiles_to_mol(smiles, verbose=False) -> Chem.Mol:
    if not verbose:
        rdBase.DisableLog('rdApp.error')
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        raise ParsingError('Could not parse SMILES string')
    if not verbose:
        rdBase.EnableLog('rdApp.error')

    return mol


def sanitize_mol(mol: Chem.Mol, verbose=False) -> Chem.Mol:
    # Ensure the charges are "standardized"
    uncharger = Uncharger()
    if not verbose:
        rdBase.DisableLog('rdApp.info')
    mol = uncharger.uncharge(mol)
    if not verbose:
        rdBase.EnableLog('rdApp.info')
    return mol


def check_charges(mol: Chem.Mol) -> None:
    # Note: formal charges are allowed but they are usually on N or O
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 and atom.GetAtomicNum() != 7 and atom.GetAtomicNum() != 8:
            logging.warning('Molecule contains charged atom that is not N or O - careful!')


def check_mol(mol: Chem.Mol):
    # Check that there aren't any hydrogen atoms left in the RDKit.Mol
    no_hs = all(atom.GetAtomicNum() != 0 for atom in mol.GetAtoms())
    if not no_hs:
        raise SanityError("Cannot process molecule: hydrogen atoms couldn't be removed")

    if NumRadicalElectrons(mol) != 0:
        raise SanityError('Molecule cannot contain radicals')

    # Check that the molecule consists of only one fragment
    fragments = Chem.GetMolFrags(mol)
    if len(fragments) != 1:
        raise SanityError(f'Incorrect number of molecular fragments ({len(fragments)})')


def embed_mol(mol, seed: int, attempt_factor=10):
    # Add hydrogen atoms in order to get a sensible 3D structure
    mol = Chem.AddHs(mol)
    # RDKit default for maxAttempts: 10 x number of atoms
    Chem.EmbedMolecule(mol, randomSeed=seed, maxAttempts=attempt_factor * mol.GetNumAtoms())
    # If not a single conformation is obtained in all the attempts, raise an error
    if mol.GetNumConformers() == 0:
        raise EmbeddingError('Generation of ligand conformation failed')
    return mol


def run_mmff94_opt(mol: Chem.Mol, max_iters: int) -> Chem.Mol:
    opt_mol = copy.copy(mol)
    Chem.MMFFSanitizeMolecule(opt_mol)
    opt_result = Chem.MMFFOptimizeMolecule(opt_mol, mmffVariant='MMFF94', maxIters=max_iters)
    if opt_result != 0:
        raise StructureOptimizationError('MMFF optimization of ligand failed')

    return opt_mol


def run_uff_opt(mol: Chem.Mol, max_iters: int) -> Chem.Mol:
    opt_mol = copy.copy(mol)
    opt_result = Chem.UFFOptimizeMolecule(opt_mol, maxIters=max_iters)
    if opt_result != 0:
        raise StructureOptimizationError('UFF optimization of ligand failed')

    return opt_mol


def refine_mol_with_ff(mol, max_iters=1000) -> Chem.Mol:
    if Chem.MMFFHasAllMoleculeParams(mol):
        try:
            opt_mol = run_mmff94_opt(mol, max_iters=max_iters)
        except Chem.rdchem.KekulizeException as exception:
            logging.info(f'Ligand optimization with MMFF94 failed: {exception}, trying UFF')
            opt_mol = run_uff_opt(mol, max_iters=max_iters)
    elif Chem.UFFHasAllMoleculeParams(mol):
        opt_mol = run_uff_opt(mol, max_iters=max_iters)

    else:
        raise StructureOptimizationError('Cannot optimize ligand: parameters not available')

    return opt_mol


def check_obabel_install():
    """ Check that openbabel is installed correctly / is the correct version """
    cmd_args = ["obabel", "-V"]
    cmd_return = subprocess.run(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout = cmd_return.stdout.decode('utf-8').strip()
    obabel_version_match = stdout == "Open Babel 3.1.0 -- Oct 12 2020 -- 14:17:21"
    if cmd_return.returncode != 0:
        raise DockingError(f"The test command `{' '.join(cmd_args)}` failed!")
    elif not obabel_version_match:
        raise DockingError("The obabel test command succeeded but the version doesn't seem to match.")


def convert_pdbqt_to_pdb(pdbqt_file: PathType, pdb_file: PathType, disable_bonding=False, verbose=False) -> None:
    # yapf: disable
    cmd_args = [
        'obabel',
        '-ipdbqt', pdbqt_file,
        '-opdb',
        '-O', pdb_file,
    ]
    # yapf: enable

    if disable_bonding:
        # "a" = read option
        # "b" = disable automatic bonding
        cmd_args += ['-ab']

    cmd_return = subprocess.run(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout = cmd_return.stdout.decode('utf-8')

    if verbose:
        logging.info(stdout)

    if cmd_return.returncode != 0:
        raise FormatConversionError('Conversion from PDBQT to PDB failed')


def protonate_mol(mol: Chem.Mol, verbose=False) -> Chem.Mol:
    smiles = Chem.MolToSmiles(mol)
    protonated_smiles = protonate_smiles(smiles, verbose=verbose)
    mol = Chem.MolFromSmiles(protonated_smiles)
    if not mol:
        raise ProtonationError(f'Cannot read protonated SMILES: {protonated_smiles}')

    return mol


def protonate_smiles(smiles: str, verbose=False) -> str:
    # Protonate SMILES with OpenBabel at given pH
    # cmd list format raises errors, therefore one string
    cmd = f'obabel -:"{smiles}" -ismi -ocan -p7.4'
    cmd_return = subprocess.run(cmd, capture_output=True, shell=True)
    output = cmd_return.stdout.decode('utf-8')

    if verbose:
        logging.info(output)

    if cmd_return.returncode != 0:
        raise ProtonationError('Ligand protonation with OpenBabel failed')

    return output.strip()


def convert_mol_file_to_pdbqt(mol_file: PathType, pdbqt_file: PathType, verbose=False):
    # yapf: disable
    cmd_list = [
        'obabel',
        '-imol', mol_file,
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
        raise FormatConversionError('Conversion from PDB to PDBQT failed')


def read_mol_from_pdb(pdb_file: PathType) -> Chem.Mol:
    mol = Chem.MolFromPDBFile(str(pdb_file))
    if not mol or mol.GetNumConformers() == 0:
        raise ParsingError(f'Cannot read PDB file {pdb_file}')
    return mol


def write_mol_to_mol_file(mol: Chem.Mol, mol_file: PathType):
    if mol.GetNumConformers() < 1:
        raise OutputError('For conversion to MDL MOL format a conformer is required')
    Chem.MolToMolFile(mol, filename=str(mol_file))


def check_vina_output(output_file: Path):
    # If Vina does not find any appropriate poses, the output file will be empty
    if os.stat(output_file).st_size == 0:
        raise DockingError('AutoDock Vina could not find any appropriate pose')


def assign_bond_orders(subject: Chem.Mol, ref: Chem.Mol, verbose=False) -> Chem.Mol:
    if not verbose:
        rdBase.DisableLog('rdApp.warning')
    try:
        mol = Chem.AssignBondOrdersFromTemplate(refmol=ref, mol=subject)
    except (ValueError, Chem.AtomValenceException) as exception:
        raise PoseProcessingError(f'Could not assign bond orders: {exception}')
    if not verbose:
        rdBase.EnableLog('rdApp.warning')
    return mol


def assign_stereochemistry(mol: Chem.Mol):
    Chem.AssignStereochemistryFrom3D(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True)


def verify_docked_ligand(ref: Chem.Mol, ligand: Chem.Mol):
    ref_smiles = Chem.MolToSmiles(ref)
    ligand_smiles = Chem.MolToSmiles(ligand)
    if ligand_smiles != ref_smiles:
        raise PoseProcessingError(f'Cannot recover original ligand: '
                                  f'{ref_smiles} (original) != {ligand_smiles} (docked)')


real_number_pattern = r'[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?'
score_re = re.compile(rf'REMARK VINA RESULT:\s*(?P<affinity>{real_number_pattern})')


def parse_affinities_from_output(output_file: PathType) -> List[float]:
    with open(output_file, mode='r') as f:
        content = f.read()
    return [float(match.group('affinity')) for match in score_re.finditer(content)]


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
