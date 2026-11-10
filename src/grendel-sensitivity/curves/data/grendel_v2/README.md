# GRENDEL v2 contour tables

The eight GRENDEL tables of v2 (2026-09-30): BC4, BC5, BC6–BC8 and BC10 with
the post-release fixes (CHANGELOG 2.0) and exHad decays (arXiv:2609.16104)
for every benchmark. A results tree in the layout the renderer reads:

```
python -m curves.plot paper --grendel-dir curves/data/grendel_v2 --benchmark hnl bc4 bc5 bc10 --out figures_v2
```

(the HNL panel needs the `curves` extra). The v1 tables stay in
`curves/data/grendel_v1`.

| File | SHA-256 |
|---|---|
| `hnl/sensitivity.csv` | `588e3403c6eb3689efb570d4f5022bcdfd07c375d5a5aa740f3746ff45d6ad8f` |
| `hnl/sensitivity_nsig10.csv` | `a902da163889eb368c72be3fa685585ab9406bcec339cb9f3869e0645f6aa597` |
| `bc4/sensitivity.csv` | `3e13ddf784907b4addab7c61ea70da2436d098a4508c0b6467c02e03fd528497` |
| `bc4/sensitivity_nsig10.csv` | `7a7948842cf39098feac3d93b7e514e98dbdf0bd907972004d521cebb9a5b8bf` |
| `bc5/sensitivity.csv` | `e608aebfcf089627f556b446bc5b486f034d31cbec5443bcd9658696afc06fa9` |
| `bc5/sensitivity_nsig10.csv` | `35b2d4f78bb10bfa3ae1bdacffd6b7f6f4f66ac5dfb1814dc72707df9acde0ae` |
| `bc10/sensitivity.csv` | `a4efa996d455ff0fe617476181877635f3e609fa480eed761022b04076906454` |
| `bc10/sensitivity_nsig10.csv` | `c94812b5f037006902b3ed41bbb364edc8d52da55775e9deb6043cbe6fb71f0c` |

`sensitivity.csv` is the scan output (`--thresholds 3 10`, both islands of
one Monte Carlo realisation); `sensitivity_nsig10.csv` is its `N_sig >= 10`
island, split off with `python -m grendel.io.thresholds ... --threshold 10`.
They are the v2 record: do not edit them, and point new scans at their own
results tree.

## How they were made

Code: this repository as tagged `v2` (the `grendel/` package is that of
536e586). Scan settings are the published ones (100 decay samples,
p > 100 MeV, the per-mass seeds) on the default mass grids.

Four-vectors:

* HNL: the v1 production with the B_c rows at `SIGMA_BC_PB` = 1.64e6 pb.
* BC4: the v1 pool (`--n-pool 400000 --high-pt-tilt-scale 5`, seed 42). It
  was produced with `particle` 0.26.2 (PDG 2024); the pinned `particle`
  1.0.1 carries b-hadron lifetimes 0.06–0.4 % shorter, so a production today
  has the same kinematics with weights 0.1–0.4 % lower per species.
* BC10: the v1 production.
* BC5: the `scan_bc5` defaults (a regeneration is byte-identical).

Decay templates: exHad 7bb2308 with Pythia 8.317, 20 000 per point, base
seed 1234, each bundle carrying exHad's own lifetime.

```
T=grendel_work/exhad      # exHad's venv Python, this package on PYTHONPATH
python -m grendel.models.hnl.templates.exhad --flavor Ue Umu Utau --out $T/hnl
python -m grendel.models.scalar.templates_exhad --out $T/scalar \
    --mass $(python -c "from grendel.models.scalar.production_bc5 import MASS_GRID_BC5 as G; print(*G)")
python -m grendel.models.alp.templates_exhad --out $T/bc10 \
    --mass $(python -c "from grendel.models.alp.mass_grid import ALP_MASS_GRID as G; print(*[m for m in G if m < 1.345])")
python -m grendel.models.alp.templates_exhad --out $T/bc10 --sub-request 250 --max-retries 10 --skip-existing
```

and the scans of the top-level README (step 4) with `--templates-dir
$T/<hnl|scalar|bc10>`.

The scalar set serves BC4 and BC5 (the BC5 grid contains the BC4 one; a
bundle depends only on its mass). BC10 was generated in two passes because
one 2000-event sub-request at 1.35 GeV failed exHad's three-body check four
times; 0.54 and 0.96 GeV are exHad's unsupported resonance windows and are
skipped (gaps in v1 too).

## Against v1 (`N_sig >= 3`, new / v1)

* HNL: within 10 % below ~2.5 GeV. At 3.4 / 3.6 GeV the lower edge is
  0.59 / 0.38 (U_e), 0.51 / 0.36 (U_μ), 0.30 / 0.27 (U_τ); exHad alone gives
  0.82 / 0.60, 0.78 / 0.55, 0.49 / 0.34 and the B_c rate the rest. U_e and
  U_μ reach 3.70 GeV (v1 3.62), U_τ 4.0 GeV (v1 3.685; peak N 3.0 there).
* BC4: lower edge 0.70–0.76 at 1.5–3.5 GeV, upper edge 0.23–0.35 at
  2–3.5 GeV: the exHad scalar-1809 width (arXiv:1809.01876 with the
  arXiv:1904.10447 corrections) is 1.85–2.3× the Winkler digitisation of v1
  there. 3.8 GeV is lost; at `N_sig >= 10` the island ends at 3.1 GeV (v1 3.6).
* BC10: lower edge within 7 %, upper edge 0.80–0.96 above 1.5 GeV; the
  marginal points 1.25–1.47 and 3.24–3.28 GeV are lost (96 sensitive, v1 103).
* BC5 (new; against its analytic-decay run of 2026-09-25): the BC4 width
  below 4 GeV (lower edge 0.49–0.59, upper 0.40–0.72 at 1–3 GeV); at
  6–20 GeV the lower edge is 1.21–1.34, the hadronic jets replacing the
  two-prong proxy.
