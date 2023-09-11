# dockstring

![CI Tests](https://github.com/mgarort/dockstring/workflows/Install%20conda%20env%20and%20run%20pytest./badge.svg?branch=main)
![Code Style: yapf](https://img.shields.io/badge/code%20style-yapf-orange.svg)

A Python package for easy molecular docking and docking benchmarking.
We can dock molecules in a few lines of code from just a SMILES string!
For details, see [our paper](https://pubs.acs.org/doi/full/10.1021/acs.jcim.1c01334)
and our [website](https://dockstring.github.io/):

> García-Ortegón, Miguel, et al. "DOCKSTRING: easy molecular docking yields better benchmarks for ligand design." Journal of Chemical Information and Modeling (2021).

## Installation

**Supported platforms:**
This package is primarily intended for Linux, but we have some support for Mac.
Please note that the scores from the Mac version do not always perfectly match the Linux version,
so we encourage the use of the Linux version whenever possible.

**Package versions:**

When installing dockstring, please be mindful of which package versions you install.
The dockstring dataset was created using:

- `rdkit=2021.03.3`
- `openbabel=3.1.1`
- `python=3.7.10`

If you want to reproduce the calculations of the dockstring dataset exactly
(or calculate docking scores completely consistent with the dataset)
then ideally install these versions of the packages above.
However, python 3.7 has reached end of life, so we have tested higher versions of the packages:
It appears that `python<=3.10, openbabel=3.1.1, rdkit<=2022.03` will also work.
Ultimately we just suggest being mindful of which version you install,
and test whether it matches the dataset values after installation (instructions on this below).
If in doubt, use our `environment.yml` file.
Note that if you do not care about consistency with our pre-computed dataset then any package version is ok.

**Installation instructions:**

We recommend installing with `conda` using our package on [conda-forge](https://anaconda.org/conda-forge/dockstring):
this will automatically install the correct versions of `rdkit` and `openbabel` (which currently cannot be installed with pip).
To do this, run:

```bash
conda install -c conda-forge dockstring
```

It can alternatively be installed from [PyPI](https://pypi.org/project/dockstring/) by running:

```bash
python3 -m pip install dockstring
```

However, this will *not* install the dependencies because `openbabel` currently cannot be installed with pip.

If you want to use dockstring for benchmarking, we recommend installing the latest version by cloning this repo:

1. Clone this repository.
1. Choose whether to install into an existing environment or create a new environment.
    - To install into a new environment, run:
      ```bash
      conda env create -f environment.yml
      conda activate dockstring
      ```
    - To install into an existing environment, simply install the desired versions of `openbabel` and `rdkit`.
1. Install the dockstring package with `pip` from this repository:
   ```bash
   pip install .
   ```
1. Check whether the installation was successful by running a test script.
   Running without error indicates a successful install.
   ```bash
   python tutorials/simple_example.py
   ```
1. *(optional)* Install [PyMol](https://pymol.org/) for target, search box and ligand visualization:
   ```bash
   conda install -c conda-forge pymol-open-source
   ```
1. *(optional)* Check whether your local version of dockstring matches the dockstring dataset.
   This is only necessary if you plan to mix pre-computed docking scores from the dockstring dataset
   with locally-computed scores, or if you want to compare results with the dockstring paper.

   We have created a `pytest` test which randomly docks `N` molecules from the dockstring dataset
   and checks whether they match. The value of `N` can be changed by setting the environment variable
   `num_dockstring_test_molecules`. We recommend starting with `N=50`, then progressing to `N=1000`
   to do a full test. The test can be run with the following commands:
   ```bash
   conda install -c conda-forge pytest  # only if not installed already
   num_dockstring_test_molecules=1000 python -m pytest tests/test_dataset_matching.py  # change "1000" to the number you wish to dock
   ```
   If the test passes then your local version of docktring matches the dataset exactly! 🥳
   If the test does not pass, we encourage you to look how the error rate (this will be displayed in the error messages).
   If 99%+ of scores match then it is probably ok to use dockstring in the benchmarks, but there will of course be some error
   and this should be noted.

If this method of installation does not work for you, please raise a github issue and we will try to help.

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
