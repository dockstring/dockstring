# Long-term reproducibility strategy for `dockstring`

## Context

`dockstring` wraps AutoDock Vina to produce *reproducible* docking scores: same SMILES in → same score out. In practice this property only holds when an exact stack is used (`rdkit=2021.03.3`, `openbabel=3.1.1`, `python=3.7.10`, vendored `vina_linux` 1.1.2 binary, Linux x86_64, default seed `974528263`, single CPU). Different OS / kernel / CPU vendor / RDKit version produces small but non-trivial drift versus the published dataset (`dockstring/resources/dataset/dockstring-dataset.tsv`, 260,148 SMILES × 58 targets, generated for the 2022 JCIM paper).

Two pressures threaten this:

1. **Stack erosion.** Python 3.7 is EOL; conda-forge will eventually stop serving `rdkit=2021.03.3` and `openbabel=3.1.1`; the bundled Vina binary is a 32-bit ELF from 2011 (last touched 2023).
2. **Native install brittleness.** The package emits a hard `DockstringError` at every `dock()` call if `obabel -V` doesn't `startswith("Open Babel 3.1")` (`dockstring/utils.py:267-280`). `pyproject.toml` declares **zero** install-time deps (`pyproject.toml:30`), so `pip install dockstring` produces a broken install. Mac scores already don't match Linux (README:17, 80).

The user has decided:
- **Effort budget**: multi-week refactor is OK.
- **Dataset stance**: keep v1 canonical; *some* regeneration is acceptable (e.g. a representative subset to quantify drift).
- **Apple Silicon**: low priority — Linux x86_64 stays the reference.
- **Image publishing**: yes, publish versioned images to a registry.

A purely technical "solve" of the determinism problem (eliminate CPU/kernel/version drift in the underlying libraries) is not feasible. RDKit's ETKDG conformer code has changed between releases; Vina's Monte Carlo search depends on FP transcendental ordering that varies across CPU microarchitectures; OpenBabel's protonation/Gasteiger-charge code has changed across 3.x. The realistic path is to **freeze the entire environment** as a published artifact and **manage drift explicitly** for users who can't run that artifact.

## Recommended approach

Treat the **published, versioned Docker image as the canonical `dockstring` environment**. Native conda becomes an officially-supported "compatibility mode" with documented drift. Add a CI gate that verifies bit-exact reproduction inside the canonical image on every PR. Capture runtime provenance so any score is self-describing. Publish a small re-docked subset to quantify how much native installs drift from v1 — this is both a service to users and a hedge for the day Python 3.7 finally falls off conda-forge.

This is a continuation and formalization of the work already started in `infra/` (the `dockstring.dockerfile` + `repro-environment.yml` + Oct 2025 `--build-arg UBUNTU_VERSION` change, commit `963c16f`). The change is from "Docker is a thing in `infra/`" to "the published Docker image *is* dockstring v1, and everything else is best-effort."

### Workstream 1 — Canonicalize and publish the image

Files:
- `infra/dockstring.dockerfile` — pin everything by exact identifier, not floating tag.
- `infra/repro-environment.yml` — already strict-pinned; add `vina` binary checksum verification.
- `infra/build_docker_image.sh` — bake in `bionic-20230530` as the only supported base for the published image; the modern-Ubuntu path becomes a separate "dev" image.
- New: `.dockerignore` at repo root — exclude the 50 MB TSV, `tmp-logs/`, `build/`, `.mypy_cache/`, `*.log` (build context is currently ~50 MB+, mostly the bundled dataset).
- New: `.github/workflows/publish-docker-image.yml` — on git tag, build and push to `ghcr.io/dockstring/dockstring:<tag>` and `:v1` (mutable v1 tag tracking the latest exact-reproduction build). Tag the image with the SHA256 of the resulting Vina binary as a label so users can verify provenance.

Specifics:
- Pin the Miniforge installer by exact version + SHA256 (currently `dockstring.dockerfile:18` fetches an unpinned URL that floats).
- Drop `--platform` flexibility — publish `linux/amd64` only. ARM64 Macs get one line of README pointing them at `--platform=linux/amd64` + Rosetta, with no fix attempt. The vendored Vina binary is x86_64-only anyway, so multi-arch images would be misleading.
- Image entrypoint: drop into the conda env so `docker run --rm -v $PWD:/work dockstring:v1 python my_script.py` Just Works.

### Workstream 2 — CI reproducibility gate

Files:
- `.github/workflows/run-tests.yml` — add a second job that runs inside the freshly-built canonical image and executes a subset of `tests/test_dataset_matching.py::test_random_matching_from_dataset` with `num_dockstring_test_molecules=200` (currently this test is `@pytest.mark.slow` and excluded from CI; see `tests/test_dataset_matching.py:60-101`).
- `tests/test_dataset_matching.py` — keep the slow marker for the full-dataset variant; add a non-slow `test_canonical_subset_matching` that runs a fixed deterministic 200-molecule subset (chosen with a seeded sample so it never changes) and is gated to only run when an env var indicates we're inside the canonical image.

This is the early-warning system. The instant a dependency upgrade or refactor breaks bit-exact reproduction inside the canonical image, the PR fails. Today there is no such alarm — the only tests in CI are `test_simple_docking` (a single `CCO`/`CC` against ABL1 in `tests/test_targets.py`) and unit tests; no real cross-target docking runs.

### Workstream 3 — Runtime provenance

Files:
- `dockstring/target.py:160-236` — extend `Target.dock` (or add a `dock_with_metadata` sister method) to return `(score, aux)` where `aux` includes:
  - `dockstring_version`
  - `rdkit_version`, `openbabel_version` (read from the existing `check_obabel_install` shellout)
  - `vina_binary_sha256` (computed once at module load, cached)
  - `platform.system()`, `platform.machine()`, `platform.libc_ver()` on Linux
  - `num_cpus` actually passed to Vina
  - `seed`
  - `canonical_environment: bool` — true iff a sentinel env var (set by the published image) is present and all version pins match
- New: `dockstring/provenance.py` — pure-functions module that builds the dict; no behavior change for existing callers.

Don't change the default return type of `dock()` (would break every downstream user). Add this as opt-in.

The `canonical_environment` flag is the user-facing payoff: a result with `canonical_environment=True` is comparable to the published dataset; `False` means it's "compatible mode" and the user has been told.

### Workstream 4 — Honest two-tier install documentation

Files:
- `README.md` — restructure the install section into two clearly-labeled tiers:
  1. **Canonical (paper-reproducible)**: `docker pull ghcr.io/dockstring/dockstring:v1`. Used for any work comparing to the published benchmarks.
  2. **Compatible (development / exploration)**: `conda env create -f environment.yml`. Documented to give scores that correlate with v1 but are not bit-exact.
- `README.md:17, 22-37, 80` — consolidate the scattered Mac/version warnings into one "About reproducibility" section that links to the drift study (Workstream 5).
- `dockstring/utils.py:267-280` — soften `check_obabel_install`: warn (don't error) for `Open Babel 3.x` versions other than 3.1.1; error only for `Open Babel 2.*` or unreadable output. Native install needs to keep working as obabel evolves.
- `pyproject.toml:30` — finally declare runtime deps with loose pins: `numpy`, `rdkit>=2021.03,<2025`, plus a clear comment that `openbabel` and `vina` (the binary) are not pip-installable. Include extras: `dockstring[paper]` could pin to the canonical versions (still won't work via pure pip, but signals intent).

### Workstream 5 — Drift characterization study

This is the "some regeneration" the user said is OK. Not a v2 dataset, just data on how much native installs drift from v1.

- New: `infra/drift_study/` — script that picks a stratified ~2000-SMILES subset of the dataset (≈35 SMILES × 58 targets, sampled deterministically) and re-docks them with the **current** rdkit/openbabel/python versions on Linux x86_64 (and once on ARM64 / Mac for completeness, even though it's not a target platform).
- Output: a small CSV of (smiles, target, v1_score, native_score, delta) plus summary statistics (R², Spearman ρ, MAE, fraction of pairs within ±0.5 kcal/mol).
- Publish results in `docs/drift.md` with a chart (or just a table). Re-run periodically (every 6 months or on major rdkit version change) and version the output.
- This becomes the empirical answer to "is native install good enough for my task?" — the user can read the drift study and decide.

Optional second pass: if the drift study shows an obvious cluster of pathological SMILES (e.g. ones that hit a changed RDKit code path), document them as known-divergent.

### Workstream 6 — Multi-CPU determinism gotcha

Files:
- `dockstring/target.py:160-236` — when called without `num_cpus` set, default to `1` and document that anything > 1 is non-deterministic even with a fixed seed (this is a Vina property, not a dockstring bug).
- Update `tests/test_targets.py::test_simple_docking` and the slow tests to assert this default.
- README "About reproducibility" section: one sentence on the multi-CPU caveat.

This is small but closes a real footgun. Currently nothing prevents a user from calling `dock(smiles, num_cpus=8)` and getting irreproducible scores while assuming the seed protected them.

### Workstream 7 — Plan-but-don't-execute v2

Files:
- `docs/v2-trigger-criteria.md` (new, terse) — the conditions under which a v2 dataset would be considered:
  - `python=3.7.10` no longer installable from conda-forge for a major Linux distro, OR
  - `rdkit=2021.03.3` no longer installable, OR
  - drift study shows R² < 0.95 vs v1 on the canonical subset, OR
  - a security issue is found in the vendored Vina binary.
- Document the v2 plan as: regenerate the entire dataset on the then-current versions, publish as a separate Docker image tag (`:v2`), keep `:v1` available indefinitely, and have `dockstring` ship both.

The point is to have the contingency in writing. Even if v2 never happens, the team has a clear, agreed-on trigger.

## Critical files to modify (when implementation begins)

- `infra/dockstring.dockerfile` — Workstream 1
- `infra/build_docker_image.sh` — Workstream 1
- `infra/repro-environment.yml` — Workstream 1
- `.github/workflows/run-tests.yml` — Workstream 2
- `.github/workflows/publish-docker-image.yml` (new) — Workstream 1
- `tests/test_dataset_matching.py:60-101` — Workstream 2
- `dockstring/target.py:160-236` — Workstreams 3, 6
- `dockstring/provenance.py` (new) — Workstream 3
- `dockstring/utils.py:267-280` — Workstream 4
- `pyproject.toml:25-30` — Workstream 4
- `README.md` (lines 17, 22-37, 80, plus install section) — Workstream 4
- `infra/drift_study/` (new) — Workstream 5
- `docs/drift.md` (new) — Workstream 5
- `docs/v2-trigger-criteria.md` (new) — Workstream 7
- `.dockerignore` (new, repo root) — Workstream 1
- `tmp_main.py`, `Session.vim`, `tmp-logs/`, `build_dockstring:*.log` — clean up untracked files left from Mac docker debugging (not part of the strategy, just hygiene)

## Suggested sequencing

If implemented in order, each workstream produces user value on its own:

1. **Workstream 1** (image publication) — gives users a one-line `docker pull` for paper-reproducible scores. Highest leverage, foundational for everything else.
2. **Workstream 2** (CI gate) — protects the canonical image from silent regressions. Should land right after #1 so the image is verified from day one.
3. **Workstream 4** (honest docs) — closes the gap between "what the README says" and "what actually reproduces."
4. **Workstream 6** (CPU determinism) — small change, big footgun closed.
5. **Workstream 3** (provenance) — moderate engineering, large user-trust payoff.
6. **Workstream 5** (drift study) — quantifies what "compatible mode" means.
7. **Workstream 7** (v2 trigger) — write-once, low effort, sets up the future.

## What this approach explicitly does *not* try to do

- Eliminate underlying randomness in RDKit/OpenBabel/Vina across CPU vendors and library versions. (Not feasible; the canonical image is the realistic version of "freezing" the stack.)
- Publish multi-arch images. (Vina binary is x86_64-only; ARM64 Macs use Rosetta + `--platform=linux/amd64`.)
- Replace the v1 dataset. (Treated as canonical forever; v2 is a contingency, not a roadmap item.)
- Add a Vina library binding to remove the subprocess call. (No reproducibility benefit; just risk.)

## Verification

For each workstream, end-to-end checks:

- **WS1**: `docker pull ghcr.io/dockstring/dockstring:v1 && docker run --rm dockstring:v1 python -c "from dockstring import load_target; t = load_target('ABL1'); print(t.dock('CCO'))"` returns `-2.4` (matching the existing `test_simple_docking` assertion in `tests/test_targets.py`).
- **WS2**: open a PR that bumps `rdkit` to `2024.x` in `repro-environment.yml`; CI must fail on the new `test_canonical_subset_matching` job.
- **WS3**: `dock_with_metadata("CCO")` returns `(score, aux)` where `aux["canonical_environment"]` is `True` inside the published image and `False` outside it.
- **WS4**: native conda install on macOS no longer hard-errors when obabel ≠ 3.1.1; instead emits a warning.
- **WS5**: `python infra/drift_study/run.py` produces `docs/drift.md` with R² > 0.99 (expected) and a delta histogram.
- **WS6**: `dock("CCO", num_cpus=8)` either produces a warning ("non-deterministic with multiple CPUs") or, if you prefer, raises unless explicitly opted into.
