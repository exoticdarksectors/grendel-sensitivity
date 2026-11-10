# GRENDEL v2: reproduction of the sensitivity figures

This repository holds the GRENDEL analysis code, grendel-sensitivity v2.5.0 (`src/grendel-sensitivity`, with its own
README), and rebuilds with it every number behind the sensitivity figures of the GRENDEL v2 paper, from scratch and
on a generic x86_64 Linux machine:

> M. Citron, F. L. Redi, R. Schmitz, *Proposing GRENDEL: GalleRy ExperimeNt for Decays of Exotic LLPs at the LHC*,
> [arXiv:2609.00152v2](https://arxiv.org/abs/2609.00152), Sec. IV "Sensitivity".

| Figure | Benchmark | Published table (in `src/grendel-sensitivity/curves/data/`) | Pipeline |
|---|---|---|---|
| Fig. 4 `bc4_grendel_paper.png` | BC4, dark scalar, B -> K S | `grendel_v2_5/bc4` | `pipelines/scalar` |
| Fig. 5 `bc5_grendel_paper.png` | BC5, dark scalar with the Higgs portal (B -> X_s SS, h -> SS) | `grendel_v2_5/bc5` | `pipelines/scalar` (+ `pipelines/senscalc` for the quartic tables) |
| Fig. 6 `hnlimits_grendel_paper.png` | BC6-BC8, heavy neutral leptons (U_e, U_mu, U_tau) | `grendel_v2_3/hnl` | `pipelines/hnl` |
| Fig. 7 `bc10_grendel_paper.png` | BC10, fermion-coupled ALP | `grendel_v2_4/bc10` (+ the FASER2/SHiP curves `bc10/projections`) | `pipelines/alp`, `pipelines/senscalc` |

Fig. 3 (the h -> SS panels `hss_panel_*`) is not covered: its published curves come from code that was never made
public, and the public code differs from them by 8-32 %.

Every pipeline follows the published recipe step by step (each README cites the records of the original runs),
runs at the published statistics, and compares every intermediate with the outputs of the original runs: byte for
byte where the platform allows it; otherwise the size and nature of each difference is stated.

## Three ways to use it

1. **Look** (no software). The rebuilt figures are in `figures/`, the rebuilt tables and SensCalc curves in
   `results/`, the comparison of each figure with the paper in `figures/comparison.txt`, and the validation reports
   of every pipeline in `results/validation/`.
2. **Redraw** (minutes). Draw the four figures from the published tables, or from the rebuilt ones, and compare
   them with the paper pixel by pixel:

   ```
   bash env/bootstrap.sh --minimal        # the analysis environment only (a few minutes)
   source env/activate.sh
   bash figures/redraw.sh published       # published tables -> identical to the paper's PNGs
   bash figures/redraw.sh rebuilt         # tables in results/ -> identical to figures/
   ```
3. **Rebuild** (about two days on 12 cores). Every pipeline from the start; see [How to run it](#how-to-run-it).

## What is in here

```
env/                 user-space toolchain: bootstrap.sh, activate.sh, the lock files of every package version
src/grendel-sensitivity/          the analysis code v2.5.0, `git archive` of commit 8802f50 of its development
                                  repository, used unpatched
src/llpatcolliders-hnl-06e5a18/   the June 2026 HNL production code (v1), still behind some HNL channels
pipelines/common/    the shared CPU pool (slot_run.py, slots_status.sh), the reference-data helpers (ref.sh, ref.py)
pipelines/fonll/     the FONLL heavy-flavour (pT, y) grids every production reads
pipelines/scalar/    BC4 and BC5: four-vectors, Higgs (pT, y) histogram, exHad templates, scans, tables, Figs. 4-5
pipelines/hnl/       HNL: MadGraph W/Z and tau pools, meson channels, exHad templates, 723 scans, table, Fig. 6
pipelines/alp/       BC10: b -> s a coupling, four-vectors, exHad templates, 153 scans, table, Fig. 7
pipelines/senscalc/  Wolfram Engine 15.0 + SensCalc: FASER2/SHiP BC10 curves, decoded decay tables, quartic tables
figures/             the rebuilt figures, one folder per pipeline; collect.sh, redraw.sh, comparison.txt
results/             the rebuilt tables (tables/), SensCalc curves (senscalc/) and validation reports (validation/)
reference/           the paper's figures (paper_figures/), the list of reference files the validations read
                     (files/, MANIFEST.tsv) and how to get them (README.md)
tools/               reference manifest and bundle, export of the public tree
third_party/         external software checked out at the pins of src/grendel-sensitivity/third_party/PINS.json
work/, logs/         data and logs of the runs (not in git)
```

External software, fetched at fixed versions: exHad e41b07b, HNLCalc (fork) b1d36c6, the HeavyN UFO model,
FONLL 1.3.3, MadGraph5_aMC@NLO 3.6.6, Pythia 8.317, LHAPDF 6.5 with NNPDF40_nlo_as_01180; for `pipelines/senscalc`:
Wolfram Engine 15.0.0, SensCalc eb2c456 (and v1.4), FeynCalc 10.2.1. HNLimits 1.2.0 (Python) draws the existing
limits of Fig. 6. Every conda and Python package is installed at the exact version of the run (`env/*.lock*`;
`env/bootstrap.sh --unlocked` resolves the newest ones instead). Licences: `THIRD_PARTY.md`.

## How to run it

Requirements: x86_64 Linux, bash, curl, git, tar, bzip2; no root. About 100 GB of disk (toolchain 4 GB, plus 14 GB
for the Wolfram Engine and SensCalc; data about 65 GB) and 13 GB of RAM or more. SensCalc needs a (free) Wolfram ID
to activate the Wolfram Engine. The reproduction was made on Ubuntu 26.04 (glibc 2.43) with a Ryzen 5 PRO 7530U
(12 cores, 13 GB).

```
bash env/bootstrap.sh                          # toolchain, Python environments, external software (~15 min)
source env/activate.sh
bash pipelines/senscalc/00_install_wolfram.sh && bash pipelines/senscalc/01_fetch_sources.sh
wolframscript -activate                        # once, with your Wolfram ID
bash run_all.sh                                # every pipeline, in parallel, sharing the CPU pool
```

`run_all.sh` starts the five pipelines together; they share one pool of CPU slots (`pipelines/common/slot_run.py`,
cores minus one, first come first served), so the machine is never oversubscribed. Each pipeline can also be run
alone (`pipelines/<area>/run_all.sh`, see its README), and every step is resume-safe: rerunning skips finished work.
`bash pipelines/common/slots_status.sh` shows what holds the CPU. At the end `figures/collect.sh` gathers the
figures and refreshes `results/`.

Network access is needed in the first steps (`env/bootstrap.sh`; `pipelines/fonll` fetches the cteq66 PDF set,
`pipelines/scalar` installs `particle` 0.26.2, `pipelines/senscalc` the Wolfram Engine, SensCalc and FeynCalc).

**Reference data (optional).** No production step reads anything but the code and the software above. The
validation steps compare with three kinds of reference: the published tables (in `src/`), the paper's figures (in
`reference/paper_figures/`), and the outputs of the original runs: campaign records, event samples, exHad
templates, SensCalc runs. The last are not in this repository (2726 files, 36 GB); `reference/README.md` lists them
and says how to get them. Point `GRENDEL_REF` at them before `source env/activate.sh`. Without them every comparison that
needs them prints a `SKIP` line (collected in `work/skipped_validations.txt`) and the run carries on: the tables
and figures come out the same, and are still compared with the published tables and the paper's figures.

**Tested.** Besides the full run: a fresh copy of the repository was installed from the lock files on the same
machine (`env/bootstrap.sh` from nothing: 15 min, every package at its locked version), where the test suite of the
analysis code passes (576 tests) and `figures/redraw.sh published` reproduces the paper's figures pixel for pixel;
and `run_all.sh` was rerun on the finished work with `GRENDEL_REF=none` (every step ends, 24 comparisons skipped)
and with the reference tree (no comparison skipped). The PDFs carry a fixed date stamp (`SOURCE_DATE_EPOCH` in
`env/activate.sh`), so rerendering a figure gives the same bytes. A run on another machine has not been made.

## Results

All four figures (Figs. 4-7) are rebuilt and validated; the last pipeline (HNL) finished on 2026-11-08.
`figures/comparison.txt` (written by `figures/collect.sh`) holds the pixel comparison of every figure with the
paper's PNG.

**The renderer is exact.** The published tables, drawn here with the unpatched plotting code, give PNGs identical to
the paper's, pixel for pixel, for all four figures (`figures/redraw.sh published`). So every difference below comes
from the data, not from the plotting.

| Figure | Rebuilt table vs published table | Rebuilt figure vs the paper's PNG |
|---|---|---|
| Fig. 4 (BC4) | 86/86 masses. Every edge bit-identical at 81; at 2.4, 2.5, 3.1, 3.6, 4.2 GeV within 0.3 %. No verdict changes (75 masses sensitive at N >= 3, 69 at N >= 10) | 0.059 % of pixels differ, none by more than 19/255. All on the contour at 2.4-4.2 GeV |
| Fig. 5 (BC5) | 112/112 masses. Every edge bit-identical at 106; at the same 5 masses and 11 GeV within 1.1 %. No verdict changes | 0.068 % of pixels differ; 0.0002 % by more than 32/255 |
| Fig. 6 (HNL) | 723/723 rows. Coarse grid (312 rows) to floating-point precision except where the Linux exHad decays differ (edges within 0.8 %). Fine grid (411 rows, 2-4.2 GeV) at the Monte Carlo level: the fine W/Z events are new MadGraph samples (cross sections agree, pull rms 0.62); edges within 1.4 % below 3.4 GeV. 2 of 2169 verdicts change, both at the Umu reach tip (3.67, 3.675 GeV: peak N 3.07/3.02 vs 2.99). Every paper number identical (best \|U_e\|^2 1.11e-8, \|U_mu\|^2 1.33e-8, \|U_tau\|^2 1.48e-7) | 0.24 % of pixels differ (0.015 % by more than 32/255), along the contours above 1.6 GeV |
| Fig. 7 (BC10) | 153/153 masses, n_hits identical at all. Every edge bit-identical at 142; at the other 11 within 0.46 %. No verdict changes. Every paper number is identical: 117/73 sensitive, reach 3.5 GeV, deepest g_Y = 2.55e-6 at 1.24 GeV | 0.049 % of pixels from the table alone. 0.64 % with the rebuilt SensCalc curves too, where the FASER2 line moves by about one pixel |

Upstream inputs:

- **FONLL grids** (`pipelines/fonll`): bottom within 3.4e-7 per bin, charm within 3.1e-5. Given the same raw grids,
  the fragmentation chain is exact.
- **SensCalc** (`pipelines/senscalc`): the BC5 quartic tables are byte-identical. The FASER2 BC10 curve agrees
  within the 1-3 % Monte Carlo noise of the SensCalc run: it samples B -> K a without a fixed seed, and three Linux
  runs scatter by up to 3.1 %. The SHiP BC10 curve, converted from M. Ovchynnikov's file, is byte-identical; the
  SensCalc rerun of it here agrees within 0.6 % in both edges.

**Why the tables are not bit-identical.** It is the platform, not the recipe. The original runs were on macOS arm64
(and lxplus); this rebuild runs on Linux x86_64.

- **exHad / Pythia.** With the same seeds, the hadronisation differs at the floating-point level. A few decays
  then differ inside one sub-request of a template bundle. This happens in about 5 % of the bundles: 5 of 86 BC4
  masses, 6 of 112 BC5 masses, 11 of 151 BC10 masses. Exactly those masses move, at the Monte Carlo level. Putting
  the published rows back at those masses gives figures pixel-identical to the paper (step 10 of
  `pipelines/scalar`). For BC10 the reverse check was run too: our four-vector pools with the published templates
  give every published edge at all 153 masses (`pipelines/alp`, control scans).
- **libm.** glibc and the macOS libm differ in the last ulp of `cos`, `sinh` and `exp`. The four-vector CSVs therefore
  differ in the 9th printed digit of a few rows, so their hashes never match.
- **pandas.** pandas 3.0.6 on macOS re-read the tables with a parser that is not round-trip; pandas 3.1.1 here reads
  them exactly. So the `sensitivity_nsig10.csv` files differ in the 14th digit.
- **SensCalc FASER2.** The run is unseeded Monte Carlo (see above).

## What is regenerated and what is taken as given

**Regenerated here, and used for the figures:** every signal event sample of BC4, BC5, BC10 and HNL (including the
MadGraph W/Z and tau samples), every exHad/Pythia decay template, the detector-geometry ray casting, all 1074 scans,
the tables and the four figures, and the FASER2 curve of Fig. 7 (SensCalc on a freshly installed Wolfram Engine).
The SHiP curve of Fig. 7 is M. Ovchynnikov's SensCalc computation (private communication, 2026-10-16): it is
converted here from his file, byte-identical to the published curve, and without that file the published curve is
used. The SensCalc rerun here gives the same curve within 0.6 % (`pipelines/senscalc`).

**Regenerated as a check, while the production reads the copy committed in the analysis code.** The original runs
read these files, and only with them can the rebuilt rows be compared with the published ones bit for bit:

| Input | Rebuilt vs committed | What the rebuilt one would change |
|---|---|---|
| FONLL b and c (pT, y) grids, read by every benchmark (`pipelines/fonll`) | within 3.4e-7 (b) and 3.1e-5 (c) per bin | negligible |
| Higgs (pT, y) histogram of BC5, MadGraph + Pythia (`pipelines/scalar`, steps 02-04) | statistically the same, not the same numbers | measured at 1, 15, 50 GeV (step 11): rows within 3 %, N >= 10 edges within 1.7 %, the scan's own Monte Carlo precision |
| K+- spectrum of the HNL kaon channels (`pipelines/hnl`, step 01) | 0.17 % more kaons (one Pythia sample) | not measured |

**Taken as given.** The analysis code itself (grendel-sensitivity v2.5.0 and the June HNL code, run unchanged); the
physics tables committed with it (the dark-scalar branching ratios and lifetimes from the literature, cross-checked
against SensCalc: within 1e-4 except one node at 0.280 GeV, 38 % off; the ALP decay tables and BC5 quartic tables,
which regenerated from SensCalc come out byte-identical); the coupling C_bs of BC10 from GKOZ Table 1, as in the
paper; and the results of other experiments drawn in the figures (existing limits and projections), except the
SHiP and FASER2 curves of Fig. 7.

## Found along the way

- **A SensCalc bug on case-sensitive file systems.** `codes/LLP distribution/prod-pheno-bremsstrahlung-proton.nb`
  imports `phenomenology/<model>/production probabilities/...`, but the directory is called `Production probabilities`.
  On macOS (case-insensitive APFS) this works. On Linux the import returns `$Failed` without an error, so proton
  bremsstrahlung is silently dropped, and the BC10 upper edges come out 13 % (SHiP) to 17 % (FASER2) low.
  `pipelines/senscalc` adds symlinks for the lower-case name (`py/case_links.py`). The bug is in SensCalc eb2c456 (the
  version the paper used, pinned here); it is fixed upstream since maksymovchynnikov/SensCalc#4, merged 2026-11-07
  (e9a5a36).
- **HFS+-compressed reference files.** Many CSVs of the original runs, as copied to the reproduction machine, read
  as empty files. The validations therefore use the digests and sizes the original runs recorded.
- **The BC10 coupling.** The ALPINIST port gives C_bs = 3.518e-4, which is a cross-check only. The production uses
  `CBS_EFF` = 4.405e-4 from GKOZ Table 1, as in the published v2.2+ runs (`pipelines/alp/README.md`).
- **A wording point in the paper.** Sec. IV (BC5) says SHiP reaches "2--8 times smaller" mixings below 2.5 GeV. Both
  the published and the rebuilt `paper_numbers.txt` give 2.1-8.9.

## Compute

On the Ryzen 5 PRO 7530U (12 cores, 11 slots), all five pipelines started together on 2026-11-06:

- FONLL: 7 h 51 min.
- SensCalc: 1 h 45 min on 4 slots.
- BC4/BC5: 21.7 slot-hours of scans, 14.2 h of wall time in the shared pool.
- BC10: 1 h for the templates and 7.6 slot-hours of scans.
- HNL: about 320 slot-hours, mostly the fine grid of 411 MadGraph W/Z runs and 411 scans. It took 39 h of wall
  time (2026-11-06 19:05 to 2026-11-08 10:04), with all 11 slots once the other pipelines had finished.

## Licence and citation

The reproduction (everything outside `src/` and `reference/paper_figures/`) and the June HNL code in
`src/llpatcolliders-hnl-06e5a18` are under the BSD 3-Clause licence (`LICENSE`); `src/grendel-sensitivity` keeps its
own (BSD-3-Clause and LGPL-3.0-or-later). `THIRD_PARTY.md` lists every external program and its licence. If you use
this repository, please cite the paper (`CITATION.cff`) and the tools the results rest on
(`src/grendel-sensitivity/CITATION.cff`).
