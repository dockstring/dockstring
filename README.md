# dockstring

![CI Tests](https://github.com/mgarort/dockstring/workflows/Install%20conda%20env%20and%20run%20pytest./badge.svg?branch=main)
![Code Style: yapf](https://img.shields.io/badge/code%20style-yapf-orange.svg)

A Python package for easy molecular docking and docking benchmarking.
We can dock molecules in a few lines of code from just a SMILES string!
For details, see [our paper](https://pubs.acs.org/doi/full/10.1021/acs.jcim.1c01334)
and our [website](https://dockstring.github.io/):

> García-Ortegón, Miguel, et al. "DOCKSTRING: easy molecular docking yields better benchmarks for ligand design." Journal of Chemical Information and Modeling (2021).

## Installation

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

If this method of installation does not work for you, please raise a github issue and we will try to help.

### Optional

Install [PyMol](https://pymol.org/) for target, search box and ligand visualization:
```bash
conda install -c conda-forge pymol-open-source 
```

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
