# Development

This document outlines development guidelines, both for new contributors and for our future selves.

## Linting

We use [pre-commit](https://pre-commit.com/) to run various linting checks.
Install by running:

```bash
conda install -c conda-forge pre-commit
pre-commit install
```

## Testing

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

## Publishing new versions of the package

dockstring is available on PyPI
([link](https://pypi.org/project/dockstring/))
and conda-forge
([link](https://anaconda.org/conda-forge/dockstring)).
Both of these are managed by @AustinT.

The workflow to publish a new version of dockstring is:

1. Ensure that the state of the main branch is ready to be published. In particular:
  - Check that all tests pass
  - Check linting (`pre-commit run --all` should work)
  - Check that `CHANGELOG.md` contains all the latest changes
2. Choose a version number for the new version: `X.Y.Z`
3. Update `CHANGELOG.md` to account for the new version: specifically, move changes from the "unreleased" to a new heading of `X.Y.Z`. Change the set of tags which count as "unreleased". Commit this to the main branch.
4. Publish a new release on github with the title `X.Y.Z`.

A GitHub workflow will automatically publish this release to PyPI,
which will automatically trigger a PR to update the release on `conda-forge`.
A maintainer of the repo <https://github.com/conda-forge/dockstring-feedstock>
will need to approve the PR to ensure that the package is updated.
