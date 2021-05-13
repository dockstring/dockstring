# dockstring

A dockstring benchmark for ML applications. (TODO improve description)

## Installation

```bash
conda env create -f environment.yml
```

## Development

### Code format

We use yapf for code formatting. This is enforced with CI checks.

```bash
yapf --style=.style.yapf --in-place --recursive ./src
```
