# BC10 (fermion-coupled ALP): Fig. 7, the GRENDEL contours

Rebuilds `curves/data/grendel_v2_4/bc10/sensitivity{,_nsig10}.csv` (the table behind Fig. 7,
`bc10_grendel_paper.png`) from scratch: the b -> s a coupling, the ALP four-vectors, the exHad decay templates,
153 one-mass scans, the tables and the figure.

    source env/activate.sh
    setsid nohup bash pipelines/alp/run_all.sh > logs/alp/run_all.out 2>&1 < /dev/null &
    cat logs/alp/run_all.status          # one line per finished step (OK / FAIL, seconds)

Every step is resume-safe: rerunning `run_all.sh` (or a single step) skips finished work. Large data go to
`work/alp/` (~14 GB, mostly the four-vector pools), logs to `logs/alp/`, reports to `work/alp/validation/`.
Knobs: `ALP_K` (parallel one-slot template jobs of step 03, default 4), `ALP_J` (scan jobs in flight, one slot each,
default 5). All CPU work goes through `pipelines/common/slot_run.py`.

**Reference material.** No step reads the original runs to produce anything. The comparisons use the published table
(`src/grendel-sensitivity/curves/data/grendel_v2_4/bc10`) and the paper's PNG (`reference/paper_figures/`), both in
this repository, and outputs of the original runs under `$GRENDEL_REF` (`grendel/campaigns/bc10/...`: the stored
macOS template bundles, the v2.4 lxplus and v2.2 scan records with their input digests, the macOS ray-cast caches;
467 files, 1.05 GB, listed in `reference/files/alp.txt`). `GRENDEL_REF` is optional (`reference/README.md`): without
it `validate_pools.py`, `validate_templates.py`, `validate_geometry.py` and the control branch (`04_scan.sh control`,
`05_collect.py control`) print a SKIP line and exit 0, while the table, its comparison with the published table and
Fig. 7 are built as usual.

## The published recipe

`curves/data/grendel_v2_4/README.md` and the campaigns (under `$GRENDEL_REF`)
`grendel/campaigns/bc10/bc10_v24_20261019` (RESUME.md, `make_templates.sh`, `lxplus/job_{prod,bc10}.sh`,
`collect.py`) and `grendel/campaigns/bc10_gkoz_cbs_20261005` (`README`, `logs/production.log`), with
`curves/data/grendel_v2_2/README.md`:

* four-vectors = the v2.2 pool (GKOZ coupling, 145-mass default grid, one RNG stream, seed 42, run on macOS at
  da9efeb + `code.diff` = commit 73f1203) plus the 8 tip masses 3.42-3.58 GeV (no 3.50) of a separate lxplus run
  over `3.4 3.42 ... 3.58` (code c621279); `--n-pool 1200000 --high-pt-tilt-scale 5`;
* templates = exHad e41b07b + Pythia 8.317 on macOS, code c621279, `--n-templates 20000 --seed 1234 --workers 6
  --sub-request 250 --max-retries 10`, 151 masses (0.54 and 0.96 GeV are inside the eta/eta' windows);
* scans = one HTCondor job per mass (153), one worker, code c621279, thresholds 3 and 10, default 60 decay samples;
  `collect.py` concatenates the rows, then `python -m grendel.io.thresholds ... --threshold 10`.

### Code provenance

`src/grendel-sensitivity` is 8802f50 (v2.5.0). Checked with `git diff` in the reference repo:

* production (`grendel/models/alp/{production,model,mass_grid}.py`, `grendel/production/**`, the tracked FONLL
  grid, `grendel/io/vectors.py`): 73f1203 -> 8802f50 changes only the model.py docstring and
  `SIGMA_KAON_PB` (HNL kaon flux, unused here). The campaign `code.diff` equals commit 73f1203 for model.py.
* templates and scans (`grendel/models/alp`, `grendel/reco`, `grendel/scan.py`, `grendel/geometry`):
  c621279 -> 8802f50 changes only the text of `geometry/data/up542/{README.md,SOURCE.json,manifest.json}`
  (grendel-geometry release tag). Numbers are unchanged, which the bit-identical scan rows below confirm. Only
  the code hash in `run_metadata.json` differs. No older code is exported and no patch is applied.

## Steps

| Step | What | Inputs -> outputs | Runtime here |
|---|---|---|---|
| `01_cbs.sh` | b -> s a coupling | `tools/compute_cbs_alpinist.py`, `model.CBS_EFF` -> `work/alp/01_cbs/cbs.txt` | 1 s |
| `02_production.sh` | four-vectors (2 runs in parallel, 1 slot each) | tracked FONLL grid -> `work/alp/vectors_v22/` (143 x 92 MB), `vectors_tip/` (10), `llp_4vectors/` (151 links) | 6 min 48 s (wall) |
| `03_templates.sh [K] [mass...]` | exHad templates, run with exHad's venv python as K one-worker, one-slot jobs over disjoint mass subsets (`work/alp/templates_jobs/job<j>`), then `merge_manifests.py` and a final manifest pass | exHad e41b07b -> `work/alp/templates_exhad/` (151 bundles + MANIFEST.json) | ~50 s per mass per core (9 masses with 3 exHad workers: 2.5 min) |
| `04_scan.sh main\|control [J] [mass...]` | one scan per mass, one worker each | pools + templates -> `work/alp/scan_<variant>/bc10_<M>/` (row, metadata, ray-cast cache, `time.txt`) | 90-175 s per mass per core (depends on the machine load), 0.44 GB RSS |
| `05_collect.py main\|control` | tables + comparison with the published rows | -> `work/alp/tables_<variant>/bc10/`, `validation/rows_<variant>.txt` | seconds |
| `06_plot.sh [variant]` | Fig. 7 from our table and from the published one, pixel diff against the paper PNG | -> `work/alp/figures_<variant>/`, `figures_published/`, `validation/figure_<variant>.txt` | seconds |
| `validate_pools.py` | byte identity of the pools via the scans' recorded input digests | -> `validation/pools.txt` | 20 s |
| `validate_templates.py` | templates vs the stored macOS bundles | -> `validation/templates.txt` | ~1 s per mass |
| `validate_geometry.py` | event-level check of the pools via the stored macOS ray-cast caches (18 masses) | -> `validation/geometry.txt` | seconds |
| `diag_ulp_sensitivity.py` | diagnostic: how many CSV rows a 1-ulp libm difference changes | scratch output | 1 min |

`run_all.sh` runs 01, 02 and validate_pools. It then runs two branches in parallel. The main branch is 03,
validate_templates, 04 main, 05 main and 06. The control branch is 04 control, 05 control and validate_geometry.
The `control` scan uses our pools with the stored macOS templates of the published run, so it isolates the
four-vectors and the platform from the templates; it needs `$GRENDEL_REF` and is skipped without it.

## Validation

**01, the coupling.** The ALPINIST port gives `cbs_t` = 3.518383e-4, the value of v1/v2, not 4.405e-4.
`model.CBS_EFF` = C_bs f_GKOZ / (m_b - m_s) = 1.8e-9 x 1e6 / (4.18 - 0.0934) = 4.404640e-4 comes from GKOZ Table 1
(CHANGELOG 2.2). The port is a cross-check (0.80x in amplitude). (CBS_EFF/ALPINIST)^2 = 1.567236, the x1.567236
of the v2.2 README. BR(B+ -> K+ a) at f_GKOZ = 1 PeV is 1.50e-7 at m_a = 0.3 GeV, as GKOZ Fig. 4.

**02, four-vectors.** The production log matches the stored v2.2 log line for line: all 143 `sigma_ref`
values (4 digits) and every count. The tip controls are 3.4 GeV -0.066 % and 3.5 GeV +0.044 % against v2.2, as
in `CERN_LOG.md`.

* *Bytes: not identical* (0 of 151). The published CSVs cannot be read on this disk: they are HFS+-compressed
  and appear empty (`campaigns/README.md`). The scans' recorded input digests (v2.4 lxplus records and the v2.2
  GKOZ record) differ for every mass. The digest construction is right: the eta/eta' pole rows, which have no
  inputs, reproduce their recorded digests.
* *Why.* The CSVs print 9 significant digits of quantities built from `cos`, `sin`, `cosh`, `sinh` and `exp`.
  macOS libm and lxplus (AVX-512 SVML in numpy) differ from glibc in the last ulp. `diag_ulp_sensitivity.py`
  moves 10 % of these results by one ulp, which changes 17-27 of the 1.2 M rows per file. So every file's hash
  differs, while the parsed values differ only at the 1e-9 level, in a few rows.
* *Same events.* The control scans reproduce the published rows bit for bit: at 1.24 GeV every field is
  identical, including n_hits = 11699; at 3.42 GeV (tip pool) peak_N differs by 1 ulp. See `rows_control.txt`
  for all 153 masses and `geometry.txt` for the per-event hit flags of 18 masses.

**03, templates.** The 9 control masses (0.3, 1.0, 1.24, 1.49, 2.0, 2.5, 3.0, 3.42, 3.5 GeV) were compared
(`validation/templates_control.txt`; all 151 in `templates.txt`):

* `ctau_m_u2eq1`, the reference ctau, the exHad table ctau, the seed and the decay backend are bit-identical. The
  mu mu convention factor is 16.51839836375364, as the macOS MANIFEST.
* Decay by decay, the same PDG sequence appears in 8 of 9 masses. At 3.0 GeV, 31 of 20 000 decays differ, all
  inside one sub-request. The largest daughter |delta p| is 3e-8 GeV and is typically 1e-16 GeV.
* exHad sub-request retries happen at the same places as on macOS (3.0 GeV: sub-request 61; 3.5 GeV:
  sub-request 3).
* All 88 homogeneity tests have p = 1.0: channel classes, labels, multiplicities, and |p| per species with a KS
  test on |p| rounded to 1e-6 GeV.

So the Linux templates are the macOS decays up to floating-point noise in the momenta.

All 151 bundles (`templates.txt`): the table fields are bit-identical in all 151. The decay sequence is identical
in 140. In the other 11 (2.15, 2.2, 2.25, 2.35, 2.65, 3.0, 3.24, 3.26, 3.46, 3.9, 4.2 GeV), a few decays differ
inside one or two sub-requests of 250 where the Pythia stream diverges; the largest daughter |dp| or |dE| over all
bundles is 3.2e-7 GeV. All 1404 homogeneity tests have p >= 0.05.

**04-06, rows and figure.** Control points at full statistics: at 1.24 GeV (the deepest point) the row computed
from our pools *and our templates* is bit-identical to the published row in all 25 fields (peak_N 6.465496186541956,
1/f = 1.0356126474856375e-8 GeV^-1, n_hits 11699). At 3.0 GeV, where 31 of the 20 000 decays differ, n_hits is
identical (10482) and the row moves at the Monte Carlo level: peak_N -0.12 %, lower edge +0.017 %, upper edge -0.09 %.
Our pools with the stored templates: 1.24 GeV bit-identical, 3.42 GeV bit-identical apart from 1 ulp in peak_N. Our tables are compared row by row (`validation/rows_main.txt`, `rows_control.txt`):
bit-identical rows, n_hits, the largest relative edge and peak_N differences, and the published numbers. The
published numbers are 117 sensitive masses, reach 3.5 GeV, deepest 1/f = 1.036e-8 GeV^-1 (g_Y = 2.55e-6) at
1.24 GeV, and the paper's g_Y ~ 2.5e-6 at 1.2-1.5 GeV. The published table renders a PNG pixel-identical to the
paper's `bc10_grendel_paper.png` (0 of 1 553 460 pixels differ).

**Final main table** (`work/alp/tables_main/bc10/`, `validation/rows_main.txt`, `figure_main.txt`; 2026-11-07):

* 153 of 153 rows; n_hits identical at all 153 masses; no mass changes verdict at N >= 3 or at N >= 10.
* Every edge (`invf_min/max`, `_open`, `_N10`) is bit-identical outside the 11 masses whose decay sequence differs
  (above). The only other differences are `peak_invf` at 1.2 GeV (4.5e-7 relative, the location of the maximum on
  the coupling grid), `peak_N` at 90 masses (<= 1.2e-14, summation order) and the effective-sample-size diagnostics
  (<= 6e-7). 45 rows are bit-identical in every field.
* At the 11 masses: `peak_N` x0.9969-1.0051 (largest at 4.2 GeV, where peak_N = 0.12, not sensitive), lower edge
  1/f x0.9954-1.0014, upper edge x0.9991-1.0039; the largest edge shifts are at 3.46 GeV, the reach tip
  (peak_N = 3.002), where the band is a sliver and its edges are ill-conditioned.
* Paper numbers, identical to the published table: 117 sensitive masses at N >= 3 and 73 at N >= 10, reach 3.5 GeV,
  deepest 1/f = 1.036e-8 GeV^-1 (g_Y = 2.55e-6) at 1.24 GeV, g_Y lower edge 2.55e-6-3.15e-6 at 1.2-1.5 GeV and
  5.24e-6-2.01e-5 below 1 GeV.
* Fig. 7 from our table (tracked SHiP/FASER2 curves): 759 of 1 553 460 pixels differ from the paper PNG (4.9e-4),
  2 by more than 32 levels (max 79), all on the contour between 2.2 and 3.5 GeV. With the rebuilt SensCalc curves
  as well (`figures/collect.sh`): 0.64 %, almost all along the FASER2 line, which moves by about one pixel (the
  1-3 % Monte Carlo noise of the SensCalc FASER2 run, `pipelines/senscalc/README.md`).
* `04_scan main` took 50 555 s (14.0 h) of wall time in the shared pool for 7.6 slot-hours of work.

**Final control table** (our pools with the stored macOS templates, `validation/rows_control.txt`; 2026-11-07 14:14):
every edge at N >= 3 and N >= 10 is bit-identical to the published table at all 153 masses, as are n_hits and every
verdict. The only differences are `peak_invf` at 1.2 GeV (4.5e-7, also in the main table, so it comes from the pools)
and the floating-point-level `peak_N` and diagnostics. So the four-vector pools reproduce the published ones exactly
as far as the scan can see, and every edge difference of the main table comes from the 11 template bundles.
Drawn as Fig. 7, the control table gives a PNG pixel-identical to the paper's (0 of 1 553 460 pixels differ,
`validation/figure_control.txt`, 2026-11-10). `validation/geometry.txt`: the ray-cast hit flags of all 1.2 M events
agree with the stored macOS caches at all 18 stored masses (entry/exit points within 3e-7 m).

`sensitivity_nsig10.csv`: the published file was written by `grendel.io.thresholds` under pandas 3.0.6 (macOS).
Its default CSV parser is not round-trip, so 573 fields of the published file differ from the values in the
published `sensitivity.csv`, by up to 7.6e-14 relative. pandas 3.1.1 here parses round-trip, so projecting the
published `sensitivity.csv` here gives exact values and a different hash. pandas 2.3.3 or 3.0.6 installed on Linux
gets all but 83 of those 573 fields (platform-dependent parser rounding). Compare this table numerically, with
rtol 1e-13.

## Status, checking and resuming

`run_all.sh` was started on 2026-11-06 at 17:55 (`setsid nohup`; status in `logs/alp/run_all.status`, step logs in
`logs/alp/`). Templates, main scans, the main table and the figure finished on 2026-11-07 at 08:58; the control scans
(our pools with the stored templates, a diagnostic that does not enter the figure), run alongside at two jobs,
finished at 14:14 (`END run_all`). Progress:

    ls work/alp/templates_jobs/*/*.npz work/alp/templates_exhad/*.npz | wc -l     # templates, of 151
    grep -c '^done' logs/alp/04_scan_main.log logs/alp/04_scan_control.log        # scan rows, of 153
    bash pipelines/common/slots_status.sh

To resume after an interruption, rerun the same command: finished pools, bundles and rows are kept.
Then read `work/alp/validation/{rows_main,rows_control,templates,geometry,figure_main}.txt`.

The 153 main scans need 7.6 slot-hours (179 s each, 0.44 GB peak memory), under an hour on 11 slots when run
alone. In the shared pool, with all five pipelines at once, each scan queued behind long HNL and BC4/BC5 jobs and
the main scans took 14 h of wall time.
