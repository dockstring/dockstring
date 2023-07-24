import copy
import logging
import os
import platform
import re
import subprocess
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Union

import pkg_resources  # type: ignore  # no type hints for this package
from rdkit import rdBase
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.Descriptors import NumRadicalElectrons
from rdkit.Chem.MolStandardize.rdMolStandardize import Uncharger

from .errors import (
    CanonicalizationError,
    DockingError,
    DockstringError,
    DockstringWarning,
    EmbeddingError,
    FormatConversionError,
    OutputError,
    ParsingError,
    PoseProcessingError,
    ProtonationError,
    SanityError,
    StructureOptimizationError,
)

PathType = Union[str, os.PathLike]


def setup_logger(level: Union[int, str] = logging.INFO, path: Optional[str] = None) -> logging.Logger:
    """
    Setup "dockstring" logger.

    :param level: log level (int or string)
    :param path: path to which log messages are written
    :return: dockstring logger
    """
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

    return logger


def is_file_empty(path) -> bool:
    """
    Check if file is empty

    :param path: path to be checked
    :return: True if path is empty, otherwise False
    """
    return os.stat(path).st_size == 0


def get_vina_filename() -> str:
    """Return system-dependent name of AutoDock Vina executable."""
    system_name = platform.system()
    if system_name == 'Linux':
        return 'vina_linux'
    if system_name == 'Darwin':
        warnings.warn(
            "Although Mac use is supported, docking scores on Mac do not always perfectly match scores from Linux. "
            "Therefore, extra care should be taken when comparing results to other platforms. "
            "In particular, the baselines in the DOCKSTRING paper were computed on Linux, "
            "so please do not directly compare your docking scores to the scores reported on the paper.",
            DockstringWarning)
        return 'vina_mac_catalina'
    else:
        raise DockstringError(f"System '{system_name}' not yet supported")


def get_resources_dir() -> Path:
    """Directory of resources (including targets and binaries)."""
    path = Path(pkg_resources.resource_filename(__package__, 'resources'))
    if not path.is_dir():
        raise DockstringError("'resources' directory not found")
    return path


def get_targets_dir() -> Path:
    """Directory containing default targets."""
    path = get_resources_dir() / 'targets'
    if not path.is_dir():
        raise DockstringError("'targets' directory not found")
    return path


def get_bin_dir() -> Path:
    """Directory containing system-dependent AutoDock Vina executables."""
    path = get_resources_dir() / 'bin'
    if not path.is_dir():
        raise DockstringError("'bin' directory not found")
    return path


def get_dataset_path() -> Path:
    """Path to dockstring dataset (not downloaded by default)."""
    return get_resources_dir() / 'dataset' / "dockstring-dataset.tsv"


def get_vina_path() -> Path:
    """Path pointing to system-dependent AutoDock Vina executable."""
    path = get_bin_dir() / get_vina_filename()
    if not path.is_file():
        raise DockstringError('AutoDock Vina executable not found')
    return path


def canonicalize_smiles(smiles: str) -> str:
    """Return new canonicalized SMILES string."""
    try:
        return Chem.CanonSmiles(smiles, useChiral=True)
    except Exception as e:
        raise CanonicalizationError(f'Cannot canonicalize SMILES: {e}')


def smiles_to_mol(smiles: str, verbose=False) -> Chem.Mol:
    """
    Attempt to convert SMILES string to RDKit Mol object

    :param smiles: SMILES string to be converted
    :param verbose: log error messages
    :return: RDKit Mol object
    """
    if not verbose:
        rdBase.DisableLog('rdApp.error')
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        raise ParsingError('Could not parse SMILES string')
    if not verbose:
        rdBase.EnableLog('rdApp.error')

    return mol


def sanitize_mol(mol: Chem.Mol, verbose=False) -> Chem.Mol:
    """
    Return new copy of molecules "standardized" charges.

    :param mol: RDKit Mol object to be standardized
    :param verbose: log error messages
    :return: new RDKit Mol object
    """
    uncharger = Uncharger()
    if not verbose:
        rdBase.DisableLog('rdApp.info')
    mol = uncharger.uncharge(mol)
    if not verbose:
        rdBase.EnableLog('rdApp.info')
    return mol


def check_charges(mol: Chem.Mol) -> None:
    """Log warning if there are on formal charges on atoms other than N or O."""
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 and atom.GetAtomicNum() != 7 and atom.GetAtomicNum() != 8:
            logging.warning('Molecule contains charged atom that is not N or O - careful!')


def check_mol(mol: Chem.Mol) -> None:
    """Check that there aren't any explicit hydrogens, radicals or multiple molecular fragments."""
    # Check that there aren't any hydrogen atoms left in the RDKit.Mol
    no_hs = all(atom.GetAtomicNum() != 1 for atom in mol.GetAtoms())
    if not no_hs:
        raise SanityError("Cannot process molecule: hydrogen atoms couldn't be removed")

    # No radicals
    if NumRadicalElectrons(mol) != 0:
        raise SanityError('Molecule cannot contain radicals')

    # Check that the molecule consists of only one fragment
    fragments = Chem.GetMolFrags(mol)
    if len(fragments) != 1:
        raise SanityError(f'Incorrect number of molecular fragments ({len(fragments)})')


def embed_mol(mol: Chem.Mol, seed: int, attempt_factor=10) -> Chem.Mol:
    """
    Attempt to generate conformation for molecule.

    :param mol: mol for which a conformation will be generated
    :param seed: seed for conformation generation
    :param attempt_factor: maximum number of attempts = <attemp_factor> * <number of atoms in molecule>
    :return: new molecule with conformation
    """
    # Add hydrogen atoms in order to get a sensible 3D structure
    mol = Chem.AddHs(mol)
    # RDKit default for maxAttempts: 10 x number of atoms
    Chem.EmbedMolecule(mol, randomSeed=seed, maxAttempts=attempt_factor * mol.GetNumAtoms())
    # If not a single conformation is obtained in all the attempts, raise an error
    if mol.GetNumConformers() == 0:
        raise EmbeddingError('Generation of ligand conformation failed')
    return mol


def run_mmff94_opt(mol: Chem.Mol, max_iters: int) -> Chem.Mol:
    """
    Optimize molecular structure with MMFF94 force field.

    :param mol: molecular structure to be optimized
    :param max_iters: maximum number of structure optimization iterations
    :return: optimized molecule
    """
    opt_mol = copy.copy(mol)
    Chem.MMFFSanitizeMolecule(opt_mol)
    opt_result = Chem.MMFFOptimizeMolecule(opt_mol, mmffVariant='MMFF94', maxIters=max_iters)
    if opt_result != 0:
        raise StructureOptimizationError('MMFF optimization of ligand failed')

    return opt_mol


def run_uff_opt(mol: Chem.Mol, max_iters: int) -> Chem.Mol:
    """
    Optimize molecular structure with UFF.

    :param mol: molecular structure to be optimized
    :param max_iters: maximum number of structure optimization iterations
    :return: optimized molecule
    """
    opt_mol = copy.copy(mol)
    opt_result = Chem.UFFOptimizeMolecule(opt_mol, maxIters=max_iters)
    if opt_result != 0:
        raise StructureOptimizationError('UFF optimization of ligand failed')

    return opt_mol


def refine_mol_with_ff(mol, max_iters=1000) -> Chem.Mol:
    """
    Optimize molecular structure. Try MMFF94 first, use UFF as a backup.

    :param mol: molecular structure to be optimized
    :param max_iters: maximum number of structure optimization iterations
    :return: optimized molecule
    """
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


def check_obabel_install() -> None:
    """Check that openbabel is installed correctly and has the correct version"""
    cmd_args = ['obabel', '-V']
    cmd_return = subprocess.run(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout = cmd_return.stdout.decode('utf-8').strip()

    if cmd_return.returncode != 0:
        raise DockingError(f"The test command `{' '.join(cmd_args)}` failed!")

    # Example: Open Babel 3.1.0 -- Oct 12 2020 -- 14:17:21
    expected_version = 'Open Babel 3.1'
    if not stdout.startswith(expected_version):
        raise DockingError("The obabel test command succeeded but the version doesn't seem to match. " +
                           expected_version + ' required.')


def convert_pdbqt_to_pdb(pdbqt_file: PathType, pdb_file: PathType, disable_bonding=False) -> None:
    """
    Convert a PDBQT file to a PDB file with Open Babel.

    :param pdbqt_file: path to the PDBQT input file
    :param pdb_file: path to the PDB output file
    :param disable_bonding: disable automatic bonding with Open Babel
    """
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
    logging.debug(stdout)

    if cmd_return.returncode != 0:
        raise FormatConversionError('Conversion from PDBQT to PDB failed')


def protonate_mol(mol: Chem.Mol, pH: float) -> Chem.Mol:
    """
    Protonate molecule at given pH

    :param mol: molecule to be protonated
    :param pH: pH at which the molecule should be protonated
    :return: protonated molecule
    """
    smiles = Chem.MolToSmiles(mol)
    protonated_smiles = protonate_smiles(smiles, pH=pH)
    mol = Chem.MolFromSmiles(protonated_smiles)
    if not mol:
        raise ProtonationError(f'Cannot read protonated SMILES: {protonated_smiles}')

    return mol


def protonate_smiles(smiles: str, pH: float) -> str:
    """
    Protonate SMILES string with OpenBabel at given pH

    :param smiles: SMILES string of molecule to be protonated
    :param pH: pH at which the molecule should be protonated
    :return: SMILES string of protonated structure
    """

    # cmd list format raises errors, therefore one string
    cmd = f'obabel -:"{smiles}" -ismi -ocan -p{pH}'
    cmd_return = subprocess.run(cmd, capture_output=True, shell=True)
    output = cmd_return.stdout.decode('utf-8')
    logging.debug(output)

    if cmd_return.returncode != 0:
        raise ProtonationError('Ligand protonation with OpenBabel failed')

    return output.strip()


def convert_mol_file_to_pdbqt(mol_file: PathType, pdbqt_file: PathType) -> None:
    """
    Convert MOL file to PDBQT file

    :param mol_file: path to MOL input file
    :param pdbqt_file: path to PDBQT output file
    """
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
    logging.debug(output)

    if cmd_return.returncode != 0 or is_file_empty(pdbqt_file):
        raise FormatConversionError('Conversion from MOL to PDBQT failed')


def read_mol_from_pdb(pdb_file: PathType) -> Chem.Mol:
    """
    Read RDKit Mol object from PDB file

    :param pdb_file: path to PDB input file
    :return: RDKit Mol object
    """
    mol = Chem.MolFromPDBFile(str(pdb_file))
    if not mol or mol.GetNumConformers() == 0:
        raise ParsingError(f'Cannot read PDB file {pdb_file}')
    return mol


def write_mol_to_mol_file(mol: Chem.Mol, mol_file: PathType) -> None:
    """
    Write RDKit Mol object to mol file

    :param mol: RDKit Mol object
    :param mol_file: Mol output path
    """
    if mol.GetNumConformers() < 1:
        raise OutputError('For conversion to MDL MOL format a conformer is required')
    Chem.MolToMolFile(mol, filename=str(mol_file))


def check_vina_output(output_file: Path) -> None:
    """
    CHeck if AutoDock Vina Output looks OK.

    :param output_file: path to AutoDock Vina output file
    """
    # If Vina does not find any appropriate poses, the output file will be empty
    if is_file_empty(output_file):
        raise DockingError('AutoDock Vina could not find any appropriate pose')


def assign_bond_orders(subject: Chem.Mol, ref: Chem.Mol, verbose=False) -> Chem.Mol:
    """
    Assign bond orders to <subject> Mol based on <ref> Mol.

    :param subject: molecules the bond orders of which are to be determined
    :param ref: reference molecule
    :param verbose: log error messages
    :return: new molecule with assigned bond orders
    """
    if not verbose:
        rdBase.DisableLog('rdApp.warning')
    try:
        mol = Chem.AssignBondOrdersFromTemplate(refmol=ref, mol=subject)
    except (ValueError, Chem.AtomValenceException) as exception:
        raise PoseProcessingError(f'Could not assign bond orders: {exception}')
    if not verbose:
        rdBase.EnableLog('rdApp.warning')
    return mol


def assign_stereochemistry(mol: Chem.Mol) -> None:
    """
    Assign stereochemistry to molecule.

    :param mol: RDKit molecule object
    """
    Chem.AssignStereochemistryFrom3D(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True)


def verify_docked_ligand(ref: Chem.Mol, subject: Chem.Mol) -> None:
    """
    Ensure that the reference and subject molecules are the same.

    :param ref: reference RDKit molecule
    :param subject: subject RDKit molecule
    """
    ref_smiles = Chem.MolToSmiles(ref)
    ligand_smiles = Chem.MolToSmiles(subject)
    if ligand_smiles != ref_smiles:
        raise PoseProcessingError(f'Cannot recover original ligand: '
                                  f'{ref_smiles} (original) != {ligand_smiles} (docked)')


real_number_pattern = r'[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?'
score_re = re.compile(rf'REMARK VINA RESULT:\s*(?P<affinity>{real_number_pattern})')


def parse_affinities_from_output(output_file: PathType) -> List[float]:
    """
    Get binding free energies (in kcal/mol) from AutoDock Vina output file.

    :param output_file: path to AutoDock Vina output file
    :return: list of binding free energies
    """
    with open(output_file, mode='r') as f:
        content = f.read()
    return [float(match.group('affinity')) for match in score_re.finditer(content)]


conf_re = re.compile(rf'^(?P<key>\w+)\s*=\s*(?P<value>{real_number_pattern})\s*\n$')


def parse_search_box_conf(conf_file: PathType) -> Dict[str, float]:
    """
    Parse search box configuration (size and position) from configuration file.

    :param conf_file: path to configuration file
    :return: dictionary containing necessary search box information
    """
    d = {}
    with open(conf_file, mode='r') as f:
        for line in f.readlines():
            match = conf_re.match(line)
            if match:
                d[match.group('key')] = float(match.group('value'))

        assert len(d) == 6
        return d
