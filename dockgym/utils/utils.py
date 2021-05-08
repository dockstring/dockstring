import os
import re
import subprocess
from typing import List, Union

from rdkit.Chem import AllChem as Chem

PathType = Union[str, os.PathLike]


class DockingError(Exception):
    """Raised when Target.dock fails at any step"""
    pass


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
        print(stdout)

    if cmd_return.returncode != 0:
        raise DockingError('Conversion from PDBQT to PDB failed')


def convert_pdb_to_pdbqt(pdf_file: PathType, pdbqt_file: PathType, verbose=False):
    # yapf: disable
    cmd_list = [
        'obabel',
        '-ipdb', pdf_file,
        '-opdbqt',
        '-O', pdbqt_file,
        '--partialcharge', 'gasteiger'
    ]
    # yapf: enable
    cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = cmd_return.stdout.decode('utf-8')

    if verbose:
        print(output)

    if cmd_return.returncode != 0:
        raise DockingError('Conversion from PDB to PDBQT failed')


def read_pdb_to_mol(pdb_file: os.PathLike):
    return Chem.MolFromPDBFile(str(pdb_file))


real_number_re = r'[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?'
score_re = re.compile(rf'REMARK VINA RESULT:\s*(?P<score>{real_number_re})')


def parse_scores_from_pdb(pdb_file: os.PathLike) -> List[float]:
    content = open(pdb_file, mode='r').read()
    return [float(match.group('score')) for match in score_re.finditer(content)]
