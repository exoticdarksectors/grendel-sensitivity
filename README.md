# grendel-sensitivity

The code that produced the GRENDEL sensitivity curves of arXiv:2609.00152v1
for the PBC benchmarks BC4 (dark scalar), BC6–BC8 (heavy neutral leptons)
and BC10 (fermiophilic ALP).

## Install

```
python -m venv .venv && . .venv/bin/activate
pip install -e ".[test]"             # Python >= 3.11
python -m pytest tests -q            # the chain end to end at tiny statistics
```

External programs are not included. `python third_party/fetch.py --list`
lists them and `--all` fetches them at the revisions in
`third_party/PINS.json`. BC4 and BC10 need only the FONLL grids in
`grendel/production/data/fonll/central/`; the HNL production needs HNLCalc
(`pip install -e ".[production]"`), MadGraph and LHAPDF; the HNL and BC10
decay templates need PyROOT with Pythia 8, and FairShip for the HNL.

## Reproduce the curves

Inputs and outputs go under `GRENDEL_WORK_DIR` (default `./grendel_work`)
unless the directories are given on the command line; every command has
`--help`.

```
# production (BC4 samples its four-vectors inside the scan)
python -m grendel.models.alp.production
python -m grendel.models.hnl.production.run_all

# decay templates (BC4 uses analytic two-body decays)
python -m grendel.models.alp.templates_pythia
python -m grendel.models.hnl.templates.fairship

# scans
python -m grendel.models.scalar.scan --thresholds 3 10 --out RESULTS/bc4
python -m grendel.models.alp.scan    --thresholds 3 10 --out RESULTS/bc10
python -m grendel.models.hnl.scan    --flavor Ue Umu Utau --thresholds 3 10 --out RESULTS/hnl
python -m grendel.io.thresholds RESULTS/<model>/sensitivity.csv --threshold 10
```

The FONLL grids are regenerated with FONLL 1.3.3 and the patches in
`grendel/production/fonll_grids/patches`
(`python -m grendel.production.fonll_grids.install --fonll DIR --build`,
then `python -m grendel.production.fonll_grids.generate`), and the kaon
spectrum with `python -m grendel.models.hnl.production.make_kaon_spectrum`.

## Citing and licence

Cite arXiv:2609.00152 and the software and data in `THIRD_PARTY.md`.
BSD-3-Clause (`LICENSE`), except the FairShip-derived files listed in
`THIRD_PARTY.md`, which are LGPL-3.0-or-later (`LICENSES/`).
