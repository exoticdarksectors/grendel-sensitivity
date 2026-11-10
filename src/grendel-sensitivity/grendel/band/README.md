# Uncertainty bands

The bands drawn as diagnostics around the exclusion islands are not
statistical confidence intervals. They are one-source-at-a-time
variations of the inputs, each re-run through the same scan as the
central curve, combined by rules that are stated here and tested in
`tests/unit`.

## Shared primitives

- `fonll_variations.py` — reads the FONLL campaign manifest
  (`variation_manifest.json` from `grendel.production.fonll_grids`),
  verifies the hash of every grid it hands out, asserts the expected counts
  (BC4 and BC10; the HNL campaign checks its bottom and charm grids with the
  same `verified_sha` and `assert_counts`) and names the variations:
  `central`, six `scale_muR*_muF*`, `pdf_0001..0100`, `mb_dn`/`mb_up`.
- `envelope.py` — the **single-source envelope** of BC4 and BC10: per
  axis the extrema of the variations (the PDF axis by its 16th/84th
  replica percentiles), the outermost endpoint over axes, numerical
  controls excluded, rebased onto the high-statistics central curve in dex.
- `quadrature.py` — the **HNL ribbon**: per-axis deviations added in
  quadrature, with a per-boundary topology veto where a variation loses the
  island or opens or closes the edge, and an explicit suppression list for points whose
  ribbon is a post-hoc policy rather than a measurement.
- `campaign_store.py` — the bookkeeping every campaign shares: config
  hashes, stable seeds, code hashes against the producer's commit,
  completion markers, compaction of finished runs, locks.

## Per model

| | variations | campaign | combine |
|---|---|---|---|
| BC4 `scalar/` | FONLL axes, the width-scheme decay alternate, two same-physics repeats (`variations.py`) | `campaign.py run/status/collect --grid-dir --scratch-dir --central-curve`; `exhad_variation.py` adds the exHad decay model with its own templates | `combine.py` → `bc4_single_source_variation_envelope.csv` |
| BC10 `alp/` | FONLL axes, `decay_gg` (gluon surrogate), `cbs` (the B → K a amplitude), two repeats (`variations.py`) | `campaign.py run/collect`, production + scan per variation as subprocesses with completion markers | `combine.py` → `bc10_single_source_variation_envelope.csv` |
| HNL `hnl/` | FONLL scale/PDF members applied coherently to bottom and charm, masses per quark (`campaign.py`) | `campaign.py --grid-dir --campaign-dir`, reusing the FONLL-independent channels | `python -m grendel.band.quadrature --registry … --out RESULTS/hnl/band/hnl_band_fonll.csv` (default `hnl_band_fonll.csv` beside the registry; `adapter.py` recombines from the raw table); `kaon_band.py` (`kaon_desc_band.csv` from the `Kmeson` runs at `d_esc` = 1, 1.5, 3 m); `decay_model_band.py` (the width-seam nuisance through lifetime and composition), `width_band.py` (its `delta(m)` table from FairShip), `bc_nuisance.py` (the B_c normalisation), `decay_model_band_exhad.py` (exHad as a third decay-model member) |

The envelopes are written next to the central curve under
`RESULTS/<model>/band/`, where `python -m curves.plot talk --with-envelope`
and `python -m curves.plot diagnostic` read them, and
`scripts/check_curves_against_bands.py` checks that a curve produced with
a different decay model sits inside its band.

## Requirements

The campaigns need the FONLL variation grids (`grendel.production.fonll_grids`,
218 grids, days of CPU), the production inputs of each model, and — for the
sizes used — the accelerated ray backend (`pip install -e ".[raycast]"`).
Every completed BC4 variation records the code state it was produced with,
and its collector refuses runs whose recorded state does not match the
named commit; the BC10 collector requires one common state across its
runs; the HNL campaign records the revision only.
