# HNL, BC6-BC8: Fig. 6

Rebuilds `curves/data/grendel_v2_3/hnl` (`sensitivity.csv`, `sensitivity_nsig10.csv`; SHA-256 `58ceddf8...`,
`2f2dc810...`), the table behind Fig. 6 (`hnlimits_grendel_paper.png`) of arXiv:2609.00152v2, from scratch on this
machine: production, MadGraph, exHad templates, scans, table, figure. Every intermediate is validated against the
stored campaign outputs (the optional reference tree `$GRENDEL_REF`), the table and figure against the published
ones (in this repository).

```
source env/activate.sh
setsid nohup bash pipelines/hnl/run_all.sh > logs/hnl/run_all.out 2>&1 &     # resume-safe; status: logs/hnl/run_all.status
```

Large data under `work/hnl/`, logs under `logs/hnl/`, validation under `work/hnl/validation/` (`summary.md`,
CSVs). Every step runs alone after `source env/activate.sh` and skips the work already done.

### Running it

`run_all.sh` needs only `env/bootstrap.sh` (toolchain, MadGraph 3.6.6, Pythia 8.317, exHad, HNLimits) and the
checked-in export `src/llpatcolliders-hnl-06e5a18`. No production step reads reference material. Step 10 compares
with the published table (in `src/`) and, where the optional reference tree `$GRENDEL_REF` is available (the
campaign outputs: `reference/README.md`), with the stored four-vectors, cross sections, template manifests and lxplus
job rows; without it those sections are reported as skipped (`GRENDEL_REF=none source env/activate.sh` runs without
it). Step 11 compares with the paper's PNG in `reference/paper_figures/`. After steps 00-03 it runs five branches in parallel, every job through the shared
slot pool (`pipelines/common/slot_run.py`): A tau pool -> coarse channels -> coarse W/Z, B templates, C fine W/Z,
D coarse scans, E fine scans. A scan point waits, holding no slot, until its template and four-vectors exist; when
the branches that produce them have ended (markers in `logs/hnl/.ended/`), a point whose inputs are still missing
fails instead of waiting for ever, and step 09 then stops with the list of missing points ("rerun run_all.sh").
Concurrent jobs per branch: `HNL_JOBS_WZ_COARSE` (2), `HNL_JOBS_WZ_FINE` (4), `HNL_JOBS_SCAN_COARSE` (2),
`HNL_JOBS_SCAN_FINE` (2). A fine scan job holds four 1M pools and peaks at 3.2 GB RSS (`time -v`, Umu 3 GeV;
~2.5 GB for most of the ray casting), so it first passes `lib/memgate.py`: at most `HNL_FINE_SCAN_MAX` (2) at once,
and only when `MemAvailable` covers its peak (3.4 GB + 1 GB margin) plus the growth still to come of the fine jobs
already running. A second driver can be started next to a running one with `lib/run_points.py ... --reverse`
(the per-point locks keep the two apart; they meet in the middle).

Measured wall time per job on this laptop (12 cores, shared with the other pipelines, one slot each):

| job | count | wall time each | slot-hours |
|---|---|---|---|
| coarse W/Z, 100k events | 312 | ~5.3 min | ~28 |
| fine W/Z, 1M events | 411 | ~25 min | ~170 |
| exHad templates, 20k decays | 723 | 0.3-7 min (mean ~1.5 min) | ~18 |
| coarse scan (~0.8-1.1M candidates) | 312 | ~2-2.5 min | ~12 |
| fine scan (4M candidates; meson pools regenerated in the job, 40 s) | 411 | ~13 min | ~90 |

i.e. ~320 slot-hours, dominated by the fine grid.

## What the published table is

`curves/data/grendel_v2_3/README.md`, campaign `grendel/campaigns/v2_3_20261016` (`collect.py`, `lxplus/`,
`hnl_fine/`) and the campaigns it reads. 723 rows = 241 masses x 3 flavours, two kinds:

| rows | four-vectors | templates | scan |
|---|---|---|---|
| **coarse**: 104 masses per flavour (repo `MASS_GRID` at 0.2-1.975 and 4.4-10 GeV; `lxplus/points_hnl.txt`) | `asbuilt_20261014/hnl/llp_4vectors/*/combined` (`lxplus/upload.sh`, `hnl_files.txt`), a mosaic of two productions (below) | `v2_baseline_20260929/templates/hnl` (= `hnl/hnl_exhad_20260917`, exHad 7bb2308) | `lxplus/job_v23.sh`, c621279 |
| **fine**: 137 masses per flavour, 2.0-4.2 GeV (`hnl/hnl_fine_condor_20261014/masses.txt`) | made inside each HTCondor job (`hnl_fine_condor_20261014/condor/job.sh`) at 10x statistics; only the SHA-256 of each combined file and the channel weight sums are stored here (the files are on EOS) | 28 masses from the v2 set, 109 generated on the Mac (`templates.sh`, `templates_exhad_new/`) | `v2_3_20261016/hnl_fine/condor/job_scan.sh`, c621279 |

The coarse four-vectors are, channel by channel ("the other channels of the v2 production"):

| channels | published production | here |
|---|---|---|
| D (charm), B_c, K (kaon) | `asbuilt_20261014/hnl/production.sh` (code 461b7e5, HNLCalc b1d36c6, Mac conda env with `particle` 0.26.2): `run_all --channels charm kaon bc --n-pool 100000 --seed 42` over the 123-mass `MASS_GRID`, then `--channels kaon --n-pool 300000 --seed 42` | step 05a |
| B mesons, Lambda_b, induced tau, prompt tau (stage 2) | the June 2026 v1 production `signal-models/hnl/llpatcolliders_HNL_PR1/hnl/tmp/runs/central_newgrids_20260623` (`production.log`): llpatcolliders 06e5a18 `hnl/run_all.py`, HNLCalc e292cef9, 116-mass grid, 100k, seed 42, all channels in one call | step 02 (code export), step 05b |
| prompt-tau pool (stage 1) | MadGraph 3.6.6, `sm`, W -> tau nu + DY tau tau, 100k events, seed offset r=21, K=1.3 (June `tmp/cache/madgraph/tau_pool`, 162 959 taus) | step 04 |
| W/Z -> l N | June v1 (`wz_fix.log`): MG5 3.6.6 + HeavyN, 100k per point, `--nb_core=1`, r=21 at every mass (the process directory was regenerated at every mass), K=1.3 | step 06 |

The fine four-vectors: `run_all --flavor F --masses M --channels bottom bc baryon --n-pool 1000000 --seed
7000/7100/7200` (B, B_c, Lambda_b), 1M MadGraph W/Z events with seed offset r=27 (Ue, Utau), r=33 (Umu), r=34 (Umu
3.845) (every `results/raw/*/generate_events_run_*.log`), the closed channels as empty files,
`combine_for_point(strict)`.

### Code versions

* `src/grendel-sensitivity` (8802f50) is used unpatched. `git diff 461b7e5 8802f50` (production) and `git diff
  c621279 8802f50` (scans) in `grendel/` touch only the two MadGraph helpers that failed on Linux (fixed in
  c621279 the way the fine campaign worked around them in `condor/compile_mg5.py`, `run_wz.py`) and documentation /
  geometry metadata; the HNL production, the scan, the geometry data and the template generator are unchanged.
* `src/llpatcolliders-hnl-06e5a18` (step 02): `git archive 06e5a18 hnl` of llpatcolliders (public at
  https://github.com/leoredi/llpatcolliders, branch `hnl-production-mesons-ew`; the export made from there is identical
  to the checked-in one), the last HNL commit before the June production (the run directory has no commit record; the outputs below confirm it). Needed because its outputs are not reproduced by the current code:
  HNLCalc e292cef9 has G_F = 1.166378e-5 (b1d36c6: 1.1663788e-5, weights x(1+1.37e-6), and the changed accept/reject
  boundary desynchronises the B-meson random stream from about 1.1 GeV on), and its SIGMA_BC_PB = 0.9e6 still feeds
  the B_c -> tau nu source of the induced-tau pool (62 vs 113 taus of 100k; the v2 B_c fix rescaled only the B_c
  channel rows). With the export, these channels come out as stored (validation 1).
* `work/hnl/venv_prod` (step 00): the production runs in a venv with the Mac production versions (numpy 2.4.6,
  scipy 1.17.1, pandas 3.0.3, `particle` 0.26.2, numba 0.65.1, sympy 1.14.0, mpmath 1.4.1; `logs/
  mac_hnl_conda_freeze.txt` of the fine campaign, which the lxplus `venv_prod` copied). Every package, transitive
  ones included, is locked in `pipelines/hnl/venv_prod.lock.txt` (installed with `pip --no-deps`). `particle` matters: HNLCalc
  takes parent lifetimes from it. The scans run in `env/venv` (numpy 2.4.6; the published scans used numpy 2.5.3,
  which has no Python 3.11 wheel).
* MadGraph 3.6.6 (`third_party`), compiled with the conda gfortran 15 (published: Homebrew gfortran on macOS for the
  June runs, the LCG gfortran on lxplus for the fine grid); HeavyN converted to Python 3 in place by step 03, once.
* exHad e41b07b (the pin) instead of 7bb2308: the two differ only in the ALP model (every HNL golden digest of
  exHad's tests is unchanged between them). Templates are made without `--reference-templates` (the FairShip ctau
  stored as a diagnostic, not read by the scan).

## Steps

| step | what | inputs -> outputs | reproduces |
|---|---|---|---|
| `00_env.sh` | production venv | PyPI, `venv_prod.lock.txt` -> `work/hnl/venv_prod` | Mac `hnl` conda env |
| `01_kaon_spectrum.sh` | rebuild the K+- SoftQCD spectrum (Pythia 8.317, 60k events, seed 42) and compare with the committed `kaon_softqcd_spectrum.npz` | -> `work/hnl/kaon_spectrum/` | `make_kaon_spectrum.py` |
| `02_export_v1_code.sh` | the June v1 production code: checks the checked-in export; without it, exports 06e5a18 from the public fork | GitHub -> `src/llpatcolliders-hnl-06e5a18` | |
| `03_mg5_procdirs.py` | MadGraph process directories (W/Z e, mu, tau; tau pool), generated and compiled once | -> `work/hnl/mg5/master/` | `ensure_process_dir` |
| `04_tau_pool.py` | prompt-tau pool, 100k, r=21 | -> `work/hnl/coarse/tau_pool.csv` | June stage 1 |
| `05_coarse_channels.sh` | coarse non-MadGraph channels (a) as-built, (b) June v1 | -> `work/hnl/coarse/{asbuilt,v1}/llp_4vectors` | see above |
| `06_wz.sh` | W/Z per point (`lib/run_wz_point.py`, driver `lib/run_points.py`): coarse 100k r=21, fine 1M r=27/33/34 (gz) | -> `work/hnl/{coarse,fine}/wz/llp_4vectors` | `run_single_point` |
| `07_templates.sh` | exHad templates, 20k, seed 1234, 723 points | -> `work/hnl/templates/` | `templates.exhad` |
| `08_scans.sh` | one scan job per point (`lib/scan_point.py`; fine jobs regenerate B/B_c/Lambda_b in seconds and combine; `lib/memgate.py`) | -> `work/hnl/results/{coarse,fine}/<fl>_<m>/` | `job_v23.sh`, `job_scan.sh` |
| `09_collect.py` | the table + N>=10 split | -> `work/hnl/table/hnl/` | `v2_3_20261016/collect.py` |
| `10_validate.py` | validation 1-5 (5: the table row by row against the published one, `lib/compare_table.py`) | -> `work/hnl/validation/` | |
| `11_plot.sh` | Fig. 6 from the rebuilt and the published table; pixel comparison with the Overleaf PNG (`reference/paper_figures`, `lib/compare_png.py`) | -> `work/hnl/figures/`, `figures/hnl/`, `work/hnl/validation/figure_*` | `curves.plot paper --benchmark hnl` |

## Validation (work/hnl/validation, `10_validate.py`; control runs in work/hnl/test)

Byte identity is out of reach on Linux wherever libm enters (last printed digit of a few rows, Pythia
hadronisation); everything else agrees exactly, so the comparisons are row by row where they can be.

| what | result |
|---|---|
| K+- spectrum rebuilt (step 01) | not identical: 10 445 of 32 000 bins differ, n_K 500 492 vs 499 629 (+0.17 %, one 60k-event Pythia sample, Linux vs macOS); sigma_inel identical. The production reads the committed file, as published |
| as-built D, B_c, K channels (step 05a) | per file: identical or 1-6 rows of 100k differing in the last printed digit (relative <= 1e-8); weights identical. Numbers per channel: `coarse_channels.csv` |
| June B, Lambda_b, induced tau (step 05b) | Ue: 8 499 813 / 8 499 917 B rows, 7 999 899 / 8 000 000 Lambda_b rows, 7 199 834 / 7 199 928 induced-tau rows text-identical to the stored blocks; the rest <= 2e-8 relative. With the current code + HNLCalc b1d36c6 instead: weights x(1+1.37e-6) and the B stream desynchronised above ~1.1 GeV; induced tau 100 000 instead of 99 999 rows (B_c rate) |
| prompt-tau pool (step 04) | sigma 5.139e4 +- 16.63 pb, 162 959 taus: as the June pool (log of `tmp/cache/madgraph/tau_pool`). Stage 2 rows: weights identical, 66 % text-identical, 98 % within 1e-3, mean energy within 6e-5 |
| coarse W/Z (step 06) | Umu 1.6 GeV, 100k, r=21: the June events. Weight sum identical (33 650.6911 pb), 66 651 / 100 000 rows text-identical, the others differ in E and p_z only (median 1.6e-8 relative: the last digits of the boost) |
| fine B, B_c, Lambda_b (scan jobs) | Umu 3 GeV, 1M pools: rows and weight sums identical to all printed digits with `production_summary.csv` (168 379.26935389556, 88 650.13400000002, 18 625.728700000003) |
| fine W/Z | Umu 3 GeV, 1M, r=33: sigma 2.585e4 +- 2.656 pb (lxplus 2.585e4 +- 2.734), weight sum 33 603.6246 vs 33 603.9366 (-9e-6); --nb_core 1 and 4 give identical files. Not the lxplus events: lxplus (LCG LHAPDF 6.5.5, gcc 14) differs from macOS (its 5k tests: +-0.16 %), and this build reproduces the macOS runs |
| templates (step 07) | ctau(U^2=1) identical. 1.0, 1.6 GeV: all 20 000 decays have the published channel and daughters, momenta equal to <= 3e-14. 2.15 GeV (Pythia strings): 19 905 / 20 000 decays with the same daughters; channel, multiplicity, momentum distributions consistent (chi2, KS p ~ 1). `lib/compare_templates.py` |
| scan (step 08) | the stored inputs of Umu 1.6 GeV (asbuilt combined file, v2 template) give the published row: every edge bit for bit, peak N to 1 ULP (summation order) |
| Fig. 6 (step 11) | the published table rendered here is pixel-identical to `current/figures/hnlimits_grendel_paper.png` (1584 x 1232; file bytes differ in PNG metadata) |

## Final results (full run, 2026-11-06 19:05 to 2026-11-08 10:04)

`work/hnl/validation/summary.md` (step 10), `figure_*.txt` (step 11). 723 of 723 rows, all scans at the published
statistics.

**Production.**
- *Coarse channels:* every channel's weight sum agrees to 1e-7 or better at all 312 points. Nearly every
  four-vector row is text-identical; the rest differ only in the last printed digit:
  - D: 15 399 668 of 15 399 866 rows;
  - B_c: 26 199 771 of 26 200 000;
  - B: 24 899 440 of 24 899 755;
  - Lambda_b: 22 299 755 of 22 300 000;
  - induced tau: 21 099 491 of 21 099 789;
  - K: 9 801 211 of 9 900 000 (the others within 9.4e-5).
- *Prompt tau and coarse W/Z:* the weight sums agree to 1e-6 or better. 66-69 % of the rows are text-identical;
  the others differ in the last digits of the boost.
- *Fine grid:* the B, B_c and Lambda_b row counts are identical at 411/411 points. The weight sums agree within
  4.4e-11; `production_summary.csv` prints B_c with 3-4 decimals only.
- *Fine W/Z, MadGraph:* the cross sections agree with lxplus at all 411 points (ratio median 1.00000, pull rms 0.62).
  The weight sums agree within 0.99995-1.00006. The events themselves are independent samples (see the fine W/Z
  row above).

**Templates.** 696 of the 723 compared with the published MANIFEST records. ctau(U^2=1) is identical at 696/696. The
channel fractions agree (no p-value below 0.01). Mean multiplicity is within 0.9978-1.0022 and the visible fraction
within 1e-4 rms.

**Table** (`work/hnl/table/hnl/`, against `grendel_v2_3/hnl`):

| block | edges (u2_min, u2_max, N >= 3 and N >= 10): median / max relative difference | n_hits identical |
|---|---|---|
| coarse (312 rows) | ~1e-12 / 0.81 % (Ue), 0.27 % (Umu), 0.24 % (Utau) | 312/312 |
| fine (411 rows) | ~1e-3 / within 1.4 % for the N >= 3 edges below 3.4 GeV; up to 21 % only where a band closes: the N >= 3 tip at 3.4-4.2 GeV and the last N >= 10 masses (3.0-3.2 GeV) | 2/411 (within 1.3 %) |

- *Coarse block:* reproduced to floating-point precision, except where the Linux exHad decays differ.
- *Fine block:* the rows move at the Monte Carlo level. The reason is that the fine W/Z four-vectors are new
  MadGraph samples, statistically equivalent to the published ones but not the same events.
- *Verdicts:* 2 of 2169 change, both at N >= 3, both in the fine Umu block at the reach tip: 3.67 GeV (peak N 3.07,
  published 2.99) and 3.675 GeV (3.02, published 2.99). The Umu reach moves from 3.665 to 3.675 GeV. No verdict
  changes at N >= 10.
- *Paper numbers, identical to the published table:*
  - best |U_e|^2 = 1.11e-8 at 0.38 GeV;
  - best |U_mu|^2 = 1.33e-8 at 1.575 GeV;
  - best |U_tau|^2 = 1.48e-7 at 3.05 GeV;
  - reach 3.66 GeV (Ue) and 3.855 GeV (Utau).

**Fig. 6.**
- *The published table rendered here:* pixel-identical to the Overleaf PNG.
- *The rebuilt table:* 4 718 of 1 951 488 pixels differ (0.24 %), 301 of them by more than 32 levels. They lie along
  the contour lines, mostly above 1.6 GeV, where the fine grid starts, and at a few coarse masses.
  `work/hnl/validation/figure_rebuilt_diff.png` marks them.
- *Copies:* `figures/hnl/hnlimits_grendel_paper.{png,pdf}`.

**Run record.** The full run took 39 h of wall time. Until 2026-11-07 09:00 it shared the pool with the other four
pipelines, and after that it had all 11 slots.

| step | finished |
|---|---|
| templates | 2026-11-07 09:40 |
| coarse W/Z | 2026-11-07 20:02 |
| coarse scans | 2026-11-07 20:10 |
| fine W/Z | 2026-11-08 07:59 |
| fine scans | 2026-11-08 10:00 |

Step 10 then failed on a pandas name clash: `x.xs` is `DataFrame.xs`, not the column `xs`. After the fix it was
rerun on 2026-11-10, in 4 min.
