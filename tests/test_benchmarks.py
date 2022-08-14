"""Simple tests for the benchmarks."""

import functools
import math
from typing import Dict

import pytest

import dockstring.benchmarks


@pytest.fixture
def original_benchmarks():
    return dockstring.benchmarks.original.get_benchmark_functions()


# For convenience, define a function for "good enough" float comparison
close_enough = functools.partial(math.isclose, rel_tol=1e-4)


class TestOriginalBenchmarks:
    """Tests for the original benchmarks."""
    def test_smoke(self, original_benchmarks: dict):
        """Can I load the benchmark functions?"""
        assert set(original_benchmarks.keys()) == {"F2", "promiscuous_PPAR", "selective_JAK2"}

    @pytest.mark.slow
    @pytest.mark.parametrize("name,top_smiles,expected_obj_parts,expected_obj", [
        ("F2", "CC1=CC2(c3cccc(C(=O)NC4=NC(=O)CC(F)=N4)c3)CC(=C1)C=CC2=O", dict(F2=-12.7, QED=0.85955), -11.2955),
        ("promiscuous_PPAR", "O=C1C=CC2=CC=C3C=C(Nc4cccc(C5C=c6cc7c(cc6=C5)CC=C7)c4)N=C3C12",
         dict(PPARA=-13.3, PPARD=-13.8, PPARG=-13.3, QED=0.829705), -11.59705),
        (
            "selective_JAK2",
            "CC1=NCC(C(N)=O)=NN1c1cc(-c2c(C)cc(C)cc2C)ccc1F",
            dict(JAK2=-12.2, LCK=-8.3, QED=0.9182718),
            -11.182718,
        ),
    ])
    def test_benchmark_values(self, original_benchmarks: dict, name: str, top_smiles: str,
                              expected_obj_parts: Dict[str, float], expected_obj: float):
        """
        Test the values of all the benchmark functions using top-scoring SMILES from the paper.
        """

        benchmark_fn = original_benchmarks[name]

        # First, check that the score formula seems correct
        # (i.e. if the parts of the score are correct, do they form the correct whole)
        assert close_enough(expected_obj, benchmark_fn.aggregation_function(**expected_obj_parts))

        # Second, check end-to end whether the SMILES produces the score expected
        actual_obj, actual_obj_parts = benchmark_fn(top_smiles)
        assert close_enough(actual_obj, expected_obj)  # overall value matches
        # Ensure parts are the same up to a certain tolerance
        assert set(actual_obj_parts.keys()) == set(expected_obj_parts.keys())
        for key in actual_obj_parts.keys():
            assert close_enough(actual_obj_parts[key], expected_obj_parts[key])

    @pytest.mark.slow
    @pytest.mark.parametrize(
        "name,nan_smiles",
        [
            ("F2", "S=C1N(C(O)C(=C(C2C(=O)N(C(=S)N(C2=O)C)C)C=3C=CN=CC3)C(=O)N1C)C"),
            ("promiscuous_PPAR", "OC1=C(C(N[C@H]2[C@H](N=C(C3=C(O)C=CC=C3)C)CCCC2)=C)C=CC=C1"),  # not NaN for first one
            ("promiscuous_PPAR", "O=C1C2(CN3C4(N(C2)CC1(C3)C(C)C)C5=C(NC4=O)C=CC=C5)C(C)C"),  # nan for all
            ("selective_JAK2", "S=C1N(C(O)C(=C(C2C(=O)N(C(=S)N(C2=O)C)C)C=3C=CN=CC3)C(=O)N1C)C"),  # nan for JAK2 only
            ("selective_JAK2", "O=C1C2(CN3C4(N(C2)CC1(C3)C(C)C)C5=C(NC4=O)C=CC=C5)C(C)C"),  # nan for LCK only
        ])
    def test_benchmark_nan(self, original_benchmarks: dict, name: str, nan_smiles: str):
        """
        Make sure everything is ok for SMILES with NaN docking scores.
        The NaN smiles are sourced from the dockstring dataset itself.
        """
        benchmark_fn = original_benchmarks[name]

        # Test 1: does it run ok?
        value, _ = benchmark_fn(nan_smiles)

        # Test 2: is the value actually nan?
        assert math.isnan(value)
