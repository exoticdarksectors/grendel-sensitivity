# The `grendel` package

One detector, one reconstruction, one scan driver; three models plugged
into it.

```
constants.py     L_int = 3 ab^-1, N_sig >= 3, the log10 |U|^2 grid, the CMS origin
geometry/        grendel_geometry.py (the fiducial mesh), reco_common.py (the two-track
                 vertex reconstruction and selection), raycast.py (batch ray-cast + cache)
reco/            acceptance.py (build_event_mc, scan_u2), exclusion.py (island solvers),
                 templates.py (decay-template bundles, exHad helpers), exhad_generation.py
scan.py          MassPoint, ScanConfig, ModelSpec, run_point, run_scan
cutflow.py       the selection cutflow from the same Monte Carlo as the scan
production/      constants, the FONLL grid parser and meson sampler, the decay kinematics,
                 data/fonll/central (the tracked grids), fonll_grids/ (their generator)
models/hnl       BC6-8: production channels (HNLCalc, MadGraph, Pythia), templates
                 (exHad, FairShip), HNLSpec, plots, channel breakdown
models/scalar    BC4: Winkler widths, B -> K S production, analytic or exHad decays, ScalarSpec
models/alp       BC10: SensCalc-derived widths and channels, B -> K a production, exHad /
                 Pythia / proxy templates, ALPSpec
band/            envelope.py (single-source envelope), quadrature.py (HNL ribbon), campaign_store.py,
                 fonll_variations.py; per model: variations, campaign, combine (+ HNL nuisances)
io/              paths.py (environment variables), vectors.py (four-vector CSVs), tables.py
                 (transparent .gz), atomic.py (durable writes), thresholds.py
```

## The scan contract

A model is a `ModelSpec` (`grendel/scan.py`): it names its mass points,
locates the four-vector file and the geometry cache of a point, supplies
the decay backend (a template bundle or an analytic engine), fixes the
random-number policy (`rng_for`, `select_events`, `chunks`) and finishes a
row from the yield curve (`finish`). `run_point` does the same work for
every model, in the same order: load four-vectors → ray-cast (cached) →
build the acceptance Monte Carlo per event chunk → `scan_u2` over the
coupling grid → solve the island at the primary threshold and at every
secondary one → diagnostics → the model's row. `run_scan` parallelises
over points, checkpoints every row atomically and writes
`run_metadata.json`.

The policies that fix the published numbers live in each model's
`spec.py` and are deliberately not unified: the HNL seeds each point by an
MD5 of (flavour, label, salt) and shares one generator across chunks; BC4
seeds by grid index, never chunks events, and keeps a 25 000-decay
reconstruction chunk inside its analytic backend; BC10 chunks 256 events
at a time with a per-chunk generator and weights templates by the visible
fraction and the matrix elements.

## Where inputs and outputs live

`grendel.io.paths` resolves everything through environment variables with
documented defaults under `GRENDEL_WORK_DIR` (`./grendel_work`):

| Variable | Meaning |
|---|---|
| `GRENDEL_<HNL/BC4/BC10>_VECTORS_DIR` | four-vector CSVs, one per mass point |
| `GRENDEL_<...>_TEMPLATES_DIR` | rest-frame decay templates (`.npz`) |
| `GRENDEL_<...>_ANALYSIS_DIR` | scan output |
| `GRENDEL_<...>_GEOMETRY_DIR` | ray-cast cache (default: `<analysis>/geometry_cache`) |
| `GRENDEL_RESULTS_DIR` | the tree the renderer reads (`<model>/sensitivity.csv`) |
| `GRENDEL_FONLL_BOTTOM_GRID`, `GRENDEL_FONLL_CHARM_GRID` | override the tracked central FONLL grids |
| `GRENDEL_FONLL_GRID_DIR` | the FONLL variation campaign (grids + manifest) |
| `GRENDEL_FONLL_DIR`, `GRENDEL_LHAPDF_DATA`, `GRENDEL_LHAPDF_CONFIG` | FONLL tree, LHAPDF data and config |
| `GRENDEL_THIRD_PARTY_DIR` | where `third_party/fetch.py` checks out |
| `GRENDEL_EXHAD_PYTHON` | the interpreter of exHad's virtual environment |
| `GRENDEL_MG5_EXE`, `GRENDEL_PYTHIA8_CONFIG` | MadGraph and Pythia executables |
| `GRENDEL_HNL_{CACHE,MG5_WORK}_DIR`, `GRENDEL_HNL_TAU_POOL_CSV` | HNL production intermediates |
| `GRENDEL_<BC4/BC10/HNL>_CAMPAIGN_DIR`, `GRENDEL_BC4_EXHAD_CAMPAIGN_DIR` | band campaign scratch |
| `GRENDEL_TRACK_P_CUT` | the track momentum cut (GeV; default 0.1) |

Command-line options override the environment. Nothing is written into the
repository at run time.

## Four-vector and template formats

Four-vector CSVs are headerless `weight, E, px, py, pz` rows (pb per event,
GeV), one file per mass point named `mN_<label>.csv` / `mS_<label>.csv` /
`mA_<label>.csv` with `<label>` = `1p000` for 1.000 GeV
(`grendel.io.vectors.format_mass_for_filename`). Template bundles are
`.npz` files with the flat arrays `daughter_counts`, `pdg`, `px`, `py`,
`pz`, `energy`, `mass`, `charge`, `stable` and the scalars `mass_GeV`,
`ctau_m_u2eq1` (the lifetime at unit coupling), `n_templates`, `seed`;
every generator writes the same bundle so the acceptance code reads one
format.
