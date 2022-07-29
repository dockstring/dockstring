"""Some generally useful functions."""

from dataclasses import dataclass
from typing import Callable, Dict, Tuple

from rdkit import Chem
from rdkit.Chem import QED as qed_module

import dockstring


@dataclass
class BenchmarkObjective:
    """
    General class for an objective function whose calculation requires
    1) evaluating a set of independent base functions
    2) aggregating those scores together into a single score
    """
    base_functions: Dict[str, Callable[[str], float]]
    aggregation_function: Callable[..., float]

    def _eval_base_functions(self, smiles: str) -> Dict[str, float]:
        return {name: f(smiles) for name, f in self.base_functions.items()}

    def __call__(self, smiles: str) -> Tuple[float, Dict[str, float]]:
        """Call all in 1."""
        base_fn_vals = self._eval_base_functions(smiles)
        return self.aggregation_function(**base_fn_vals), base_fn_vals


def safe_dock_function(smiles: str, target_name: str, **dock_kwargs):
    """Call dockstring and return nan if there are errors."""
    target = dockstring.load_target(target_name)
    try:
        docking_output = target.dock(smiles, **dock_kwargs)
        score = docking_output[0]
    except dockstring.DockstringError:
        score = float("nan")
    return score


def QED(smiles: str) -> float:
    """Calculates QED from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    return qed_module.qed(mol)
