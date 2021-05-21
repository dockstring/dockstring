import math
import os
from pathlib import Path

from dockstring.utils import parse_single_pdbqt

resources_dir = Path(os.path.dirname(os.path.realpath(__file__))) / 'resources'


def test_parse_count():
    with open(resources_dir / 'ligand.pdbqt', mode='r') as f:
        content = f.read()

    atoms = parse_single_pdbqt(content)

    assert len(atoms) == 22

    assert atoms[0].symbol == 'O'
    assert math.isclose(atoms[0].positions[0], -1.724)

    assert atoms[-1].symbol == 'C'
    assert math.isclose(atoms[-1].positions[-1], -0.554)
