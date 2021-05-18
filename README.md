# dockstring

![CI Tests](https://github.com/mgarort/dockstring/workflows/Install%20conda%20env%20and%20run%20pytest./badge.svg?branch=main)
![Code Style: yapf](https://img.shields.io/badge/code%20style-yapf-orange.svg)

A dockstring benchmark for ML applications. (TODO improve description)

## Installation

```bash
conda env create -f environment.yml
```

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
