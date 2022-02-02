# dockstring

![CI Tests](https://github.com/mgarort/dockstring/workflows/Install%20conda%20env%20and%20run%20pytest./badge.svg?branch=main)
![Code Style: yapf](https://img.shields.io/badge/code%20style-yapf-orange.svg)

A Python package for easy molecular docking.
For details, see our paper:

**DOCKSTRING: easy molecular docking yields better benchmarks for ligand design**<br>
Miguel García-Ortegón, Gregor N. C. Simm, Austin J. Tripp, José Miguel Hernández-Lobato, Andreas Bender, Sergio Bacallado<br>
https://arxiv.org/abs/2110.15486

Our **dataset** containing docking scores and poses of more than 260K ligands for 58 medically-relevant targets
can be downloaded [here](https://figshare.com/s/95f2fed733dec170b998?file=30562257).

## Installation

1. Create a new `conda` environment from the `environment.yml` file in this repository:
   ```bash
   conda env create -f environment.yml
   ```
1. Activate the new `conda` environment:
   ```bash
   conda activate dockstring
   ```
1. Install the dockstring package with `pip`:
   ```bash
   pip install .
   ```

### Optional

Install [PyMol](https://pymol.org/) for target, search box and ligand visualization:
```bash
conda install -c conda-forge pymol-open-source 
```

## Tutorials

- Check out this [tutorial](tutorials/1_docking_risperidone_against_DRD2.ipynb) to get started!
- This [tutorial](tutorials/tbd.ipynb) shows you how to access our dataset.

## Development

### Code format

We use yapf and flake8 for code formatting.
Run the following to check formatting:

```bash
yapf --style=.style.yapf --in-place --recursive .
flake8 --config=.flake8 .
```

We have CI set up to check this, but we _highly_ recommend setting up
[pre-commit](https://pre-commit.com/) to avoid accidentally committing bad code.
You can do so in the following way:

```bash
conda install -c conda-forge pre-commit
pre-commit install
```
