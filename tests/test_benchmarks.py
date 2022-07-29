"""Simple tests for the benchmarks."""

import pytest

import dockstring.benchmarks.original


@pytest.fixture
def original_benchmarks():
    return dockstring.benchmarks.original.get_benchmark_functions()


class TestOriginalBenchmarks:
    def test_smoke(self, original_benchmarks: dict):
        """Can I load the benchmark functions?"""
        assert set(original_benchmarks.keys()) == {"F2", "promiscuous_PPAR", "selective_JAK2"}

    def test_f2_value(self, original_benchmarks: dict):
        """Does F2 match a known SMILES?"""
        assert True  # TODO


# TODO
# match some known values
# check that NaNs are ok
