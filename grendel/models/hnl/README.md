# Heavy neutral leptons (BC6–BC8)

One HNL mixing with a single active flavour; `|U|²` is the scan variable,
`m_N` runs over `mass_grid.MASS_GRID`, and the lifetime and visible final
states come from the decay templates, not from the production rates.

## Production

Each channel writes `<vectors>/<flavor>/<channel>/mN_<label>.csv` and
`combine_channels` sums them into `<flavor>/combined/`, which the scan
reads. `production/run_all.py` runs the whole set in parallel.

| Channel | Module | Needs |
|---|---|---|
| `Dmeson`, `Bmeson`, `Bc` — meson → ℓ N (2- and 3-body) | `generate_meson_csvs` | HNLCalc; the FONLL charm/bottom grids; the B_c rate `SIGMA_BC_PB` |
| `Bbaryon` — Λ_b → Λ_c ℓ N | `generate_baryon_csvs` | HNLCalc |
| `tau` — prompt τ → N X from W and Drell–Yan | `madgraph/run_tau_production` (the τ pool), `tau_decay` | MadGraph + LHAPDF |
| `induced_tau` — τ from D_s decays | `generate_induced_tau` | HNLCalc |
| `Kmeson` — K → ℓ N | `generate_kaon_csvs` | the tracked Pythia 8.315 SoftQCD kaon spectrum (`make_kaon_spectrum` regenerates it) |
| `WZ` — W → ℓ N, Z → ν N | `madgraph/run_wz_production`, `run_wz_sharded`, `lhe_to_csv` | MadGraph, the HeavyN UFO model, LHAPDF |

Production is done at unit coupling; the scan applies `|U|²` explicitly.
`transport_control.py` is the transport cross-check of the kaon channel.
The kaon channel's escape-distance nuisance (`--d-esc`, default
`KAON_D_ESC`) is a variation of `generate_kaon_csvs`: regenerate the
`Kmeson` channel at the alternative distances, recombine, rescan, and
compare the lower edges.

## Decay templates

`templates/exhad.py` generates the exHad rest-frame decays (the model of
the current curves) through `grendel.reco.exhad_generation`, with the
width table exHad supplies giving `ctau(|U|² = 1)`; it must run under the
interpreter of exHad's virtual environment (`GRENDEL_EXHAD_PYTHON`).
`templates/fairship.py` generates the FairShip decays of the first
published curve under PyROOT + Pythia 8, using FairShip's `python/`
modules from `third_party/fetch.py FairShip` and the channel selection in
`data/fairship_decay_selection.conf`. Both write the common bundle format.

## Scan and plots

```
python -m grendel.models.hnl.scan --flavor Ue Umu Utau --thresholds 3 10
python -m grendel.models.hnl.plot_exclusion RESULTS/hnl/sensitivity.csv
python -m grendel.models.hnl.channel_breakdown        # which channel carries the lower edge
python -m grendel.models.hnl.plot_money --results-dir RESULTS/hnl   # the curve with its bands
```

`spec.py` holds the numerical policy: `seed_for(flavor, label, salt)`
(MD5), one generator shared across event chunks, exact hit selection
unless `--max-hit-events` resamples, the grid solver
`find_exclusion_band`, and the width-delta extra scans the decay-model
band uses. `data/width_band_delta.csv` is the seam-derived `delta(m)`
table (`grendel.band.hnl.width_band`), precomputed because computing it
needs FairShip under PyROOT.
