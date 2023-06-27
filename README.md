# dockstring

![CI Tests](https://github.com/mgarort/dockstring/workflows/Install%20conda%20env%20and%20run%20pytest./badge.svg?branch=main)
![Code Style: yapf](https://img.shields.io/badge/code%20style-yapf-orange.svg)

A Python package for easy molecular docking and docking benchmarking.
We can dock molecules in a few lines of code from just a SMILES string!
For details, see [our paper](https://pubs.acs.org/doi/full/10.1021/acs.jcim.1c01334)
and our [website](https://dockstring.github.io/):

> García-Ortegón, Miguel, et al. "DOCKSTRING: easy molecular docking yields better benchmarks for ligand design." Journal of Chemical Information and Modeling (2021).

## Installation

This package is primarily intended for Linux, but we have some support for Mac.
However, the scores from the Mac version does not match the Linux version,
so we encourage the use of the Linux version whenever possible.

To ensure compatibility with the dockstring dataset,
the package has very strict versioning requirements
for its main dependencies (`rdkit` and `openbabel`).
As such, we recommend you install it in the following way.

1. Clone this repository.
1. Choose whether to install into an existing environment or create a new environment.
    - To install into an existing environment, install the correct version of `openbabel` and `rdkit`:
      ```bash
      conda install -c conda-forge rdkit=2021.03.3 openbabel=3.1.1
      ```
    - To install into a new environment, run:
      ```bash
      conda env create -f environment.yml
      conda activate dockstring
      ```
1. Install the dockstring package with `pip` from this repository:
   ```bash
   pip install .
   ```
1. Check whether the installation was successful by running a test script.
   Running without error inducates a successful install.
   ```bash
   python tutorials/simple_example.py
   ```
1. *(optional)* Install [PyMol](https://pymol.org/) for target, search box and ligand visualization:
   ```bash
   conda install -c conda-forge pymol-open-source 
   ```

If this method of installation does not work for you, please raise a github issue and we will try to help.
## Matching dockstring dataset

One usage of dockstring is to run the benchmarks associated with the dockstring dataset.
Although dockstring is designed to produce the same numbers on different platforms,
in pratice we find that on some systems the docking scores returned by dockstring may differ.
If you plan to mix docking scores from the dockstring dataset with locally-computed scores,
or if you want to compare results with the dockstring paper, we *strongly* encourage you to
check whether your local version of dockstring matches the dockstring dataset.

We have created a `pytest` test which randomly dock `N` molecules from the dockstring dataset
and checks whether they match. The value of `N` can be changed by setting the environment variable
`num_dockstring_test_molecules`. We recommend starting with `N=50`, then progressing to `N=1000`
to do a full test. The test can be run with the following commands:
```bash
conda install -c conda-forge pytest  # only if not installed already
num_dockstring_test_molecules=1000 python -m pytest tests/test_dataset_matching.py  # change "1000" to the number you wish to dock
```
If the test passes then your local version of docktring matches the dataset! 🥳

If the test does not pass, we encourage you to look how the error rate (this will be displayed in the error messages).
If 99%+ of scores match then it is probably ok to use dockstring in the benchmarks, but there will of course be some error
and this should be noted.

## Tutorials

- See dockstring's basic usage [here](tutorials/1_docking_risperidone_against_DRD2.ipynb).
- Learn how to visualize docking poses [here](tutorials/2_visualizing_dataset_poses.ipynb)

See [our website](https://dockstring.github.io/) for linkks to tutorials for
our dataset and benchmarks.

## Development

We use [pre-commit](https://pre-commit.com/) to enforce code formatting and style.
Install by running:

```bash
conda install -c conda-forge pre-commit
pre-commit install
```

We use [pytest](https://docs.pytest.org) to test our code.
You can install pytest by running `conda install -c conda-forge pytest`.
Before committing, please run the following to make sure that all tests pass:

```bash
python -m pytest tests/
```

Alternatively, to skip a variety of slow tests, run:

```bash
python -m pytest -m "not slow" tests/
```

## Citation

If you use the dockstring package/dataset/benchmark in your work,
please use the following citation:

```tex
@article{garciaortegon2022dockstring,
    author = {García-Ortegón, Miguel and Simm, Gregor N. C. and Tripp, Austin J. and Hernández-Lobato, José Miguel and Bender, Andreas and Bacallado, Sergio},
    title = {DOCKSTRING: Easy Molecular Docking Yields Better Benchmarks for Ligand Design},
    journal = {Journal of Chemical Information and Modeling},
    volume = {62},
    number = {15},
    pages = {3486-3502},
    year = {2022},
    doi = {10.1021/acs.jcim.1c01334},
    URL = {https://doi.org/10.1021/acs.jcim.1c01334}
}
```
