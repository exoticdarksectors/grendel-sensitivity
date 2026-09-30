# grendel-sensitivity

The code that produces the GRENDEL sensitivity curves of arXiv:2609.00152
for the PBC benchmarks BC4 (dark scalar), BC5 (dark scalar with
BR(h → SS) = 0.01), BC6–BC8 (heavy neutral leptons) and BC10 (fermiophilic
ALP).

## Install

```
python -m venv .venv && . .venv/bin/activate
pip install -e ".[test]"             # Python >= 3.11
python -m pytest tests -q            # the chain end to end at tiny statistics
```

External programs are not included. `python third_party/fetch.py --list`
lists them and `--all` fetches them at the revisions in
`third_party/PINS.json`. BC4, BC5 and BC10 need only the FONLL grids in
`grendel/production/data/fonll/central/` (BC5 also the Higgs histogram in
`grendel/models/scalar/data/`); the HNL production needs HNLCalc
(`pip install -e ".[production]"`), MadGraph and LHAPDF; the HNL and BC10
decay templates need PyROOT with Pythia 8, and FairShip for the HNL; the
exHad templates run in exHad's own virtual environment (with Pythia 8.317).

## Reproduce the curves

Inputs and outputs go under `GRENDEL_WORK_DIR` (default `./grendel_work`)
unless the directories are given on the command line; every command has
`--help`.

The commands below give the v2 curves (tag `v2`); the v1 curves of
arXiv:2609.00152v1 come from tag `v1`.

```
# production (BC4 and BC5 sample their four-vectors inside the scan)
python -m grendel.models.alp.production
python -m grendel.models.hnl.production.run_all

# exHad decay templates: exHad's venv Python, this package on PYTHONPATH
T=grendel_work/exhad
python -m grendel.models.hnl.templates.exhad --flavor Ue Umu Utau --out $T/hnl
python -m grendel.models.scalar.templates_exhad --out $T/scalar \
    --mass $(python -c "from grendel.models.scalar.production_bc5 import MASS_GRID_BC5 as G; print(*G)")
python -m grendel.models.alp.templates_exhad --out $T/bc10 \
    --mass $(python -c "from grendel.models.alp.mass_grid import ALP_MASS_GRID as G; print(*[m for m in G if m < 1.345])")
python -m grendel.models.alp.templates_exhad --out $T/bc10 --sub-request 250 --max-retries 10 --skip-existing

# scans
python -m grendel.models.scalar.scan --n-pool 400000 --high-pt-tilt-scale 5 --templates-dir $T/scalar --thresholds 3 10 --out RESULTS/bc4
python -m grendel.models.scalar.scan_bc5 --templates-dir $T/scalar --thresholds 3 10 --out RESULTS/bc5
python -m grendel.models.alp.scan    --templates-dir $T/bc10 --thresholds 3 10 --out RESULTS/bc10
python -m grendel.models.hnl.scan    --flavor Ue Umu Utau --templates-dir $T/hnl --thresholds 3 10 --out RESULTS/hnl
python -m grendel.io.thresholds RESULTS/<model>/sensitivity.csv --threshold 10
```

The published BC4 pool was produced with `particle` 0.26.2; with the pinned
1.0.1 (PDG 2025 b-hadron lifetimes) its weights come out 0.1–0.4 % lower.

The FONLL grids are regenerated with FONLL 1.3.3 and the patches in
`grendel/production/fonll_grids/patches`
(`python -m grendel.production.fonll_grids.install --fonll DIR --build`,
then `python -m grendel.production.fonll_grids.generate`), and the kaon
spectrum with `python -m grendel.models.hnl.production.make_kaon_spectrum`.
The BC5 Higgs histogram is rebuilt with MadGraph and Pythia 8
(`grendel.models.scalar.higgs_production`, then `higgs_shower` and
`higgs_pool`).

## Citing and licence

Cite arXiv:2609.00152 and the software and data in `THIRD_PARTY.md`.
BSD-3-Clause (`LICENSE`), except the FairShip-derived files listed in
`THIRD_PARTY.md`, which are LGPL-3.0-or-later (`LICENSES/`).
