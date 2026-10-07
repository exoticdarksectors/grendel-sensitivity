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
arXiv:2609.00152v1 come from tag `v1`. Tag `v2.1` adds the per-production-mode
and rescaled-BR(h → SS) islands to the BC5 table (below); tag `v2.2` changes
the BC10 production coupling and gives the curves of arXiv:2609.00152v2
(below).

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
python -m grendel.io.variants RESULTS/bc5/sensitivity.csv --variant hSS     # --list names them all
```

### The curves of arXiv:2609.00152v2 (tag `v2.2`)

```
# BC4, BC5: exHad final states with the Winkler lifetime and e e / mu mu / tau tau / hadronic split
python scripts/make_hybrid_templates.py $T/scalar $T/scalar_hybrid
python -m grendel.models.scalar.scan --n-pool 400000 --high-pt-tilt-scale 5 --templates-dir $T/scalar_hybrid --thresholds 3 10 --out RESULTS/bc4
python -m grendel.models.scalar.scan_bc5 --templates-dir $T/scalar_hybrid --thresholds 3 10 --out RESULTS/bc5

# BC10: b -> s a coupling of arXiv:2310.03524 Table 1 (this tag's CBS_EFF, rates x1.567 over v2)
python -m grendel.models.alp.production --n-pool 1200000 --high-pt-tilt-scale 5
python -m grendel.models.alp.scan --templates-dir $T/bc10 --thresholds 3 10 --out RESULTS/bc10
```

The HNL curves are the v2 ones with 3.05, 3.10, 3.15, 3.25, 3.30, 3.35, 3.45,
3.50 and 3.55 GeV added and every mass of 2.6–3.7 GeV recomputed at 4x
production statistics (`run_all --channels bottom bc baryon wz --no-prompt-tau
--skip-combine --n-pool 400000 --wz-nevents 400000 --seed 4242`, the channels
closed at these masses written empty, then `combine_for_point`).

### BC5 by production mode

BC5 makes the scalar three ways: `mixing` (b → X_s S through the mixing
angle, the BC4 channel), `hSS` (on-shell h → SS, which only an LHC
experiment has) and `BSS` (b → X_s SS and B_s → SS through the off-shell
Higgs, the only Higgs-mediated piece a beam dump sees). Besides the nominal
island, every row of the BC5 `sensitivity.csv` carries the island of each
mode alone (`mixing_u2_min`, `hSS_peak_N`, `BSS_has_sensitivity`, ...) and,
for each `--br-hss-overlay` value (default 0.001), the islands rescaled to
that BR(h → SS): the full curve (`brhss0.001_*`) and each quartic mode
(`brhss0.001_hSS_*`, `brhss0.001_BSS_*`). The quartic yields are linear in
BR(h → SS), so all of these come from the same Monte Carlo as the nominal
curve. `python -m grendel.io.variants <sensitivity.csv> --variant hSS` writes
`sensitivity_hSS.csv` in the layout of `sensitivity.csv`, and
`grendel.io.thresholds` applies to it as to any table.

The mode of every four-vector row is the seventh column of the BC5 CSVs. CSVs
written by v2 lack it; the scan regenerates them (same seeds, same pools).

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
