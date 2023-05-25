"""Tests that the output of the dock function matches the values in the dataset."""

from __future__ import annotations

import math
import os
import random
import urllib.request

import pytest

from dockstring import load_target

# List of random target/value pairs from dockstring dataset
random_dataset_tuples = [
    ('PDE5A', 'O(CCC)C(=O)C=1C=C(OC)C(O)=CC1', -6.9),
    ('HSD11B1', 'C=1(C=CC(C[C@H](N)C(=O)N2CSCC2)=CC1)C3=C(C=C(C=C3)F)F', -9.1),
    ('PPARA', 'O(C(=O)C(C(C1=CC=CC=C1)C2=CC=CC=C2)=CC(=O)O)CC', -7.9),
    ('SRC', 'C1(=C(C(Cl)=CC(=C1)Cl)OCC(CC(CC(=O)O)O)O)C(C2=CC=C(C=C2)F)C3=CC=C(C=C3)F', -9.2),
    ('DRD2',
     'N(C([C@H](CCC(=O)O)NC([C@H](CC(=O)O)NC([C@H](C)NC([C@H](CC(=O)O)N)=O)=O)=O)=O)[C@H](CC1=CC=CC(C(P(=O)(O)O)(F)F)=C1)C(N[C@H](C(N)=O)CC(C)C)=O',
     -8.9),
    ('ESR1', 'O=C1C(NC2=CC=C(N)C=C2)=C(NC(=O)CC)C(=O)C3=C1C=CC=C3', -8.6),
    ('ADORA2A', 'IC1=CC(F)=C(NC2=C(C(OC3=CC=4NC(OC4C=C3)=O)=CC(F)=C2)C(=O)N)C=C1', -10.5),
    ('ADORA2A', 'S(C1=CC(F)=C(C=C1)CCCO)(=O)(=O)N', -6.7),
    ('PDE5A', 'C=12N(C(=O)N(C(C1N=C(N2)C3=CC(S(N4CCN(CC4)C)(=O)=O)=CC=C3OCCC)=O)C)CC(C)C', -10.1),
    ('HSP90AA1', 'O=C(N1CCN(C2CCCCC2)CC1)C=3C=4C(C(=O)N(C3)CC)=CC(OC)=C(OC)C4', -10.0),
    ('AKT2', 'N(C)(C(C1=CNC2=CC=C(C=C21)Br)C=3C=CC(=CC3)OC)C=4C=CC=CC4', -8.3),
    ('MMP13', 'BrC12CC3(CC(C1)CC(C3)C2)CC(=O)NCC4=CC=CC=C4', -9.1),
    ('EGFR', 'O(C[C@@H](O)/C=C/CC(OC)OC)CC1=CC=CC=C1', -6.9),
    ('NR3C1', 'ClC=1C(OC(=O)CC2=CC(S(=O)(=O)N3CCOCC3)=C(OC)C=C2)=CC=CC1', -8.5),
    ('MAPK14', 'S(NC=1N(N=CC1)C2=CC=CC=C2)(=O)(=O)C=3C=CC=CC3', -7.0),
    ('HSD11B1', 'C1=CC(=C(C=C1C(C(CC)NC2CCCC2)O)O)O', -8.0),
    ('HMGCR', 'S1C=C(NC(=O)NC=2C=CC(C=3C=4C(=NNC4N)N=CC3)=CC2)C=C1', -7.9),
    ('MMP13', 'S(=O)(=O)(N(C)C)C=1C=C(NC(=O)COC(=O)CCC2CCCCC2)C(=CC1)C', -8.9),
    ('MMP13', 'ClC1=CC=C(CO/N=C(/C=2SC(=CC2)CC#N)\\C)C=C1', -8.4),
    ('PGR', 'C(C=1C=C(N2CCOCC2)C=C(C(NC=3C=C(NC(C=4C=CC(OCC=5N=CC=CC5)=CC4)=O)C(=CC3)C)=O)C1)(F)(F)F', -9.0),
    ('MET', 'N1C2=C(C=CC(=C2)CN3CCN(CC3)C4=CC=C(C=C4)C)OCC1=O', -9.7),
    ('CDK2', 'N1(C2=C(C(N)=NC(=N2)OCCC3=C(Cl)C=CC=C3)N=C1)C4[C@@H]([C@H](O)[C@H](O4)CO)O', -9.0),
    ('KDR', 'ClC1=CC=C(COC=2C=CC(C=3OC(NCCCN4C=CN=C4)=C(N3)C#N)=CC2)C=C1', -9.9),
    ('MAOB', 'O1C2=C(CN(C)C)C(O)=CC=C2C(=O)C(OC3=CC=C(OC)C=C3)=C1C', -7.1),
]


@pytest.mark.slow
@pytest.mark.parametrize("target_name,smiles,value", random_dataset_tuples)
def test_value_matches(target_name: str, smiles: str, value: float):
    target = load_target(target_name)
    computed_value, _ = target.dock(smiles)
    assert isinstance(computed_value, float)
    assert math.isclose(value, computed_value)


@pytest.fixture(scope="module")
def whole_dockstring_dataset() -> list[tuple[str, str, float]]:
    """
    Loads the whole dockstring dataset as a list of (target, SMILES, score) tuples.
    """

    # Download the dataset from figshare and read it in
    file_name, _ = urllib.request.urlretrieve("https://figshare.com/ndownloader/files/35948138")
    with open(file_name) as f:
        lines = f.readlines()

    # Read the dataset line by line
    output_tuples: list[tuple[str, str, float]] = []
    header = lines[0].strip().split()
    for line in lines[1:]:
        items = line.strip().split()
        for i, docking_score in enumerate(items[2:]):
            output_tuples.append((header[i + 2], items[1], float(docking_score)))
    return output_tuples


@pytest.mark.slow
def test_random_matching_from_dataset(whole_dockstring_dataset: list[tuple[str, str, float]]) -> None:

    # Sanity check: make sure random tuples above are in the dataaset
    dataset_as_set = set(whole_dockstring_dataset)
    for t in random_dataset_tuples:
        assert t in dataset_as_set

    # Shuffle the dataset
    rng = random.Random(888)
    rng.shuffle(whole_dockstring_dataset)

    # Choose the first N molecules to dock, where N can be customized
    num_docking_scores_to_validate = int(os.environ.get("num_dockstring_test_molecules", 100))
    dataset_to_use = whole_dockstring_dataset[:num_docking_scores_to_validate]
    del whole_dockstring_dataset

    # Dock all molecules
    scores: list[float] = []
    expected_scores: list[float] = []
    for target_name, smiles, true_score in dataset_to_use:
        target = load_target(target_name)
        computed_value, _ = target.dock(smiles)
        assert isinstance(computed_value, float)
        scores.append(computed_value)
        expected_scores.append(true_score)

    # Now check that scores match
    assert len(scores) == len(expected_scores)
    scores_match = [math.isclose(s, e) for s, e in zip(scores, expected_scores)]
    if not all(scores_match):
        print("Scores do not match for the following molecules:")
        for i, (target_name, smiles, true_score) in enumerate(dataset_to_use):
            if not scores_match[i]:
                print(target_name, smiles, true_score, scores[i])
        raise AssertionError
