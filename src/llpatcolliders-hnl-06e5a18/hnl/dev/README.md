# `hnl/dev/` - internal development utilities

This directory contains exploratory or diagnostic tools that are useful while
developing the package but are not required to reproduce the production and
analysis chain.

## `debug_plots.py`

This script plots yield versus mass, channel fractions, and HNL spectra from a
completed run:

```bash
conda activate hnl
cd hnl
python dev/debug_plots.py --channel-fraction Ue
python dev/debug_plots.py --channel-fraction all
```

It reads `tmp/runs/<tag>/llp_4vectors/` and writes below
`tmp/runs/<tag>/analysis/debug_plots/`.

Production and analysis code must not import from `hnl/dev/`.

## Environment

Use the committed package environment:

```bash
conda env create -f hnl/environment.yml
conda activate hnl
python hnl/dev/debug_plots.py --help
```

For a one-shot invocation without activation:

```bash
conda run -n hnl python hnl/dev/debug_plots.py --help
```

## PR Preparation

The publishable package includes production, analysis, FairShip decay inputs,
static HNL inputs, tools, tests, and top-level entry points. The GRENDEL
geometry and reconstruction are **not** part of the `hnl/` PR -- they live in
the sibling `higgs/` package already on `origin/main` (PR #13), and the hnl PR
must not modify `higgs/`. A branch assembled from another development branch
must include all of the hnl pieces:

```bash
git fetch origin
git checkout -b pr/<topic> origin/main

git checkout hnl-production-mesons-ew -- \
    .gitignore \
    hnl/.gitignore \
    hnl/README.md \
    hnl/REMAINING_WORK.md \
    hnl/environment.yml \
    hnl/config_mass_grid.py \
    hnl/run_all.py \
    hnl/run_analysis.py \
    hnl/run_full_all.sh \
    hnl/analysis/ \
    hnl/data/ \
    hnl/production/ \
    hnl/tools/ \
    hnl/tests/ \
    hnl/vendored/

git status --short
conda activate hnl
python -P -m pytest hnl/tests/ -q
```

Do not include `hnl/tmp/`, generated plots, caches, logs, or local MadGraph
work trees. The current branch already integrates `origin/main` (PR #13) via a
merge, so this manual cherry-pick is only an alternative way to assemble a clean
hnl-only branch from scratch.

## Placement Rule

Files needed by another user to reproduce the result belong in one of:

- `hnl/production/`
- `hnl/analysis/`
- `hnl/data/`
- `hnl/tools/`
- `hnl/tests/`
- `hnl/vendored/`
- a documented top-level entry point

The GRENDEL geometry and reconstruction are not in `hnl/`; they are imported
from the shared `../higgs/` single source (PR #13 on `origin/main`).

Exploratory utilities that are not part of the reproducible chain belong in
`hnl/dev/`.
