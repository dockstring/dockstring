"""Benchmarks for original dockstring paper."""

import functools
import math
from typing import Dict

from .utils import QED as QED_function
from .utils import BenchmarkObjective, safe_dock_function


# Raw functions
def QED_penalty(qed: float) -> float:
    return 10.0 * (1.0 - qed)


def F2_score(*, F2: float, QED: float) -> float:
    return F2 + QED_penalty(QED)


def promiscuous_PPAR_score(*, PPARA: float, PPARD: float, PPARG: float, QED: float) -> float:
    # Max of a list of NaNs is not always NaN so we check this manually
    if any(math.isnan(v) for v in [PPARA, PPARD, PPARG]):
        return math.nan
    return max(PPARA, PPARD, PPARG) + QED_penalty(QED)


def selective_JAK2_score(*, JAK2: float, LCK: float, QED: float) -> float:

    # Note: there was a small error in the formula for this objective
    # in our JCIM publication. *THIS* is the correct formula,
    # which matches the numbers in all of our tables
    lck_median_score = -8.1
    return JAK2 - min(LCK - lck_median_score, 0) + QED_penalty(QED)


def get_benchmark_functions(**dock_kwargs) -> Dict[str, BenchmarkObjective]:
    """
    Returns the functions for the original benckmarks.

    dock_kwargs specifies kwargs to pass to the `target.dock` function
    (e.g. pH, num_cpus)
    """
    output: Dict[str, BenchmarkObjective] = dict()

    # F2
    output["F2"] = BenchmarkObjective(
        base_functions=dict(
            F2=functools.partial(safe_dock_function, target_name="F2", **dock_kwargs),
            QED=QED_function,
        ),
        aggregation_function=F2_score,
    )

    # Promiscuous PPAR
    ppar_funcs = {
        target_name: functools.partial(safe_dock_function, target_name=target_name, **dock_kwargs)
        for target_name in ["PPARA", "PPARD", "PPARG"]
    }
    output["promiscuous_PPAR"] = BenchmarkObjective(
        base_functions=dict(
            QED=QED_function,
            **ppar_funcs,
        ),
        aggregation_function=promiscuous_PPAR_score,
    )

    # Selective JAK2
    jak2_funcs = {
        target_name: functools.partial(safe_dock_function, target_name=target_name, **dock_kwargs)
        for target_name in ["JAK2", "LCK"]
    }
    output["selective_JAK2"] = BenchmarkObjective(
        base_functions=dict(
            QED=QED_function,
            **jak2_funcs,
        ),
        aggregation_function=selective_JAK2_score,
    )

    return output
