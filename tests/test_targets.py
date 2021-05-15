import pytest

from dockstring import list_all_target_names, load_target, DockingError


class TestLoader:
    def test_load_all_targets(self):
        names = list_all_target_names()
        assert all(isinstance(name, str) for name in names)
        assert all(load_target(name) for name in names)

    def test_wrong_target(self):
        with pytest.raises(DockingError):
            load_target('does_not_exist')
