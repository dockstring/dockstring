import os
import subprocess
import tempfile
from pathlib import Path
from typing import List

import pkg_resources
from rdkit.Chem import AllChem as Chem

from dockgym.dockgym import Target
from dockgym.utils import parse_search_box_conf, write_embedded_mol_to_pdb


def view(target: Target, mols: List[Chem.Mol] = None, search_box=True):
    """
    Start pymol and view the receptor and the search box.
    """
    commands = ['pymol', target._pdb]

    if search_box:
        pymol_view_search_box_file = pkg_resources.resource_filename(__package__,
                                                                     os.path.join('utils', 'view_search_box.py'))
        conf = parse_search_box_conf(target._conf)
        # yapf: disable
        commands += [
            pymol_view_search_box_file,
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
