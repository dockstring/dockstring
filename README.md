# dockstring

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
