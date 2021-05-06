# Exceptions


class DockingError(Exception):
    """Raised then Target.dock fails at any step"""
    pass


# Chemistry


def get_num_conf(mol):
    return sum(1 for conformer in mol.GetConformers())
