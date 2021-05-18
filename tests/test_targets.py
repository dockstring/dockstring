import math

import pytest

from dockstring import list_all_target_names, load_target, DockingError


class TestLoader:
    def test_load_all_targets(self):
        names = list_all_target_names()
        assert len(names) == 58
        assert all(isinstance(name, str) for name in names)
        assert all(load_target(name) for name in names)

    def test_wrong_target(self):
        with pytest.raises(DockingError):
            load_target('does_not_exist')


class TestDocking:
    def test_simple_docking(self):
        target = load_target('ABL1')

        smiles_1 = 'CCO'
        energy_1, _ = target.dock(smiles_1)
        assert math.isclose(energy_1, -2.4)

        smiles_2 = 'CC'
        energy_2, _ = target.dock(smiles_2)
        assert math.isclose(energy_2, -1.8)
