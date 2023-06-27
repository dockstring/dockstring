"""Tests that the output of the dock function matches the values in the dataset."""

from __future__ import annotations

import math
import os
import random

import pytest

from dockstring import load_target, dataset as dockstring_dataset


@pytest.fixture()
def whole_dockstring_dataset() -> list[tuple[str, str, float]]:
    """
    Loads the whole dockstring dataset as a list of (target, SMILES, score) tuples.
    """

    output_tuples: list[tuple[str, str, float]] = []
    dataset = dockstring_dataset.load_dataset()
    for target, score_dict in dataset.items():
        for smiles, score in score_dict.items():
            output_tuples.append((target, smiles, score))
    return output_tuples


def test_dataset_values_correct(whole_dockstring_dataset: list[tuple[str, str, float]]) -> None:
    """
    Tests that the values in the dockstring dataset match some hard-coded ones
    and that the right number of values are present.
    """
    dataset_as_set = set(whole_dockstring_dataset)
    for t in dockstring_dataset.random_dataset_tuples:
        assert t in dataset_as_set
    assert len(whole_dockstring_dataset) == 58 * 260148  # 58 targets, 260148 unique SMILES in total


@pytest.mark.slow
def test_random_matching_from_dataset(whole_dockstring_dataset: list[tuple[str, str, float]]) -> None:
    """
    Randomly docks molecules from the dataset and checks that the scores match.

    The number of molecules is controlled by the environment variable num_dockstring_test_molecules.
    If not set, it defaults to 100.
    """

    # Shuffle the dataset
    rng = random.Random(888)
    rng.shuffle(whole_dockstring_dataset)

    # Choose the first N molecules to dock, where N can be customized
    num_docking_scores_to_validate = int(os.environ.get("num_dockstring_test_molecules", 100))
    dataset_to_use = whole_dockstring_dataset[:num_docking_scores_to_validate]
    del whole_dockstring_dataset

    # Dock all molecules and check whether scores match
    observed_scores: list[float] = []
    non_matching_points = []
    for target_name, smiles, expected_score in dataset_to_use:
        target = load_target(target_name)
        observed_score, _ = target.dock(smiles)
        assert isinstance(observed_score, float)
        observed_scores.append(observed_score)

        # Matching
        if not math.isclose(observed_score, expected_score):
            non_matching_points.append((target_name, smiles, expected_score, observed_score))
    assert len(observed_scores) == len(dataset_to_use)  # sanity check

    # Test: assert that all scores match and give detailed error message if they don't
    max_deviation = max([0.] + [abs(s - e) for _, _, s, e in non_matching_points])
    error_str = (f"Scores do not match for {len(non_matching_points)}/{len(dataset_to_use)} molecules. " +
                 f"Max deviation: {max_deviation}. Scores without match:")
    for target_name, smiles, expected_score, observed_score in non_matching_points:
        error_str += f"\nTarget: {target_name}\nSMILES: {smiles}\nExpected score: {expected_score}\nObserved score: {observed_score}\n"
    assert len(non_matching_points) == 0, error_str
