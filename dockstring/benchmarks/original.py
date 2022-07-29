"""Benchmarks for original dockstring paper."""

import functools
from typing import Dict
from .utils import BenchmarkObjective, safe_dock_function, QED as QED_function


# Raw functions
def QED_penalty(qed: float) -> float:
    return 10.0 * (1.0 - qed)


def F2_score(*, F2: float, QED: float) -> float:
    return F2 + QED_penalty(QED)


def promiscuous_PPAR_score(*, PPARA: float, PPARD: float, PPARG: float, QED: float) -> float:
    return max(PPARA, PPARD, PPARG) + QED_penalty(QED)


def selective_JAK2_score(*, JAK2: float, LCK: float, QED: float) -> float:
    return JAK2 - min(LCK, -8.1) + QED_penalty(QED)


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
