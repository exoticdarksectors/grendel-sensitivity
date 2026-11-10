# BC4 and BC5 (dark scalar): Figs. 4 and 5

Rebuilds `curves/data/grendel_v2_5/{bc4,bc5}` (the tables behind Figs. 4 and 5 of arXiv:2609.00152v2) from
scratch on this machine, following the published recipe step by step, and validates every intermediate against
the stored campaign outputs.

```
source env/activate.sh
setsid nohup bash pipelines/scalar/run_all.sh > logs/scalar/run_all.out 2>&1 &     # resume-safe; status in logs/scalar/run_all.status
```

Large data under `work/scalar/`, logs under `logs/scalar/`, validation reports under `work/scalar/validation/`
(summaries quoted below). Every step can be run alone (`bash pipelines/scalar/NN_*.sh`) after
`source env/activate.sh`; each skips the work already done. `run_all.sh` runs 00, 01, 05, 06, then 07 and 08 while
02 -> 03 -> 04 run alongside, then 09, 10 and 11; a failing step stops it (see `run_all.status`), a rerun resumes.

Needs, besides `env/bootstrap.sh`: network once (00 installs `particle` 0.26.2 from PyPI).

**Reference material** (`pipelines/common/ref.sh`; the productions read none of it). The published tables
(`src/grendel-sensitivity/curves/data/grendel_v2_5`) and the paper's PNGs and PDFs (`reference/paper_figures`) are in
this repository, so the comparisons with them (07, 09, 10) always run. The outputs of the original runs (campaign
records, the stored 2026-09-25 LHE files and their shower, the stored BC5 four-vectors, the published exHad bundles)
are the optional tree `$GRENDEL_REF` (`reference/README.md`; the files this pipeline reads are
listed in `reference/files/scalar.txt`). Without it every step still runs to the end: each comparison that needs
it prints a `SKIP` line (collected in `work/skipped_validations.txt`) and the step's OK line counts them. Skipped
then: the size check of 01, the MG5 check of 02, the shower of the stored LHE in 03 (and with it (b) of 04), the
vector check of 05, the template check of 06 (so 09 compares all masses in one group and 10 draws no spliced
figures), and the `control_reftpl` scans of 07.

## The published recipe (what v2.5 is)

`curves/data/grendel_v2_5/README.md` + `grendel/campaigns/scalar/scalar_eventcalc_20261019` (RESUME.md,
`lxplus/job_sc.sh`, `collect.py`): the v2.3 runs (code c621279, geometry 2026-10-16.1, the v2.3 four-vectors of
`v2_3_20261016/lxplus/upload.sh`, one HTCondor job per mass, one worker, thresholds 3 and 10, BC5 with
`--decay-samples 100`) with the decay templates replaced by the pure exHad `scalar-1809` bundles of
`v2_baseline_20260929/templates/scalar` (no Winkler hybrid). The inputs trace back to:

| input | published origin | here |
|---|---|---|
| BC4 four-vectors | `scratch/bc4_nsig10_20260820` (Aug 2026, grendel-mono): `--n-pool 400000 --high-pt-tilt-scale 5`, seed 42, the 86 masses of `production.MASS_GRID` in one RNG stream, `particle` 0.26.2 (`grendel_v2/README.md`) | step 01 |
| BC5 four-vectors | `asbuilt_20261014/chain.sh` (code 461b7e5): `scan_bc5` defaults (n_pool 200 000, n_higgs 200 000, seed 42), the 112 masses in the interleaved order of the chain, one RNG stream; tracked FONLL grid, SensCalc quartic tables and Higgs histogram | step 05 |
| Higgs (pT, y) histogram | `grendel/models/scalar/data/higgs_pt_y_14TeV.csv.gz` (built 2026-10-14 from `asbuilt_20261014/bc5/madgraph`: MG5 LO of 2026-09-25 + Pythia 8.317 shower) | steps 02-04 |
| exHad templates | `v2_baseline_20260929/templates/scalar`: exHad 7bb2308, Pythia 8.317, 20 000 per mass, base seed 1234, sub-requests of 2000 with 3 retries (86 BC4 masses from `bc4_exhad_20260918`, 26 BC5 masses from `logs/templates_bc5.log`) | step 06 |
| scans | `scalar_eventcalc_20261019/lxplus/job_sc.sh` | steps 07-08 |
| tables, Sec. IV numbers | `collect.py`, `bc5_v2_1_20261003/results/bc5_hybrid/paper_numbers.py`, `v2_3_20261016/edge_ratios.py` (copied to `py/ref/`, unchanged apart from the default `CURVES_ROOT` of `paper_numbers.py` and the anonymised `rsync` line in the docstring of `collect.py`) | step 09 |
| figures | `python -m curves.plot paper --grendel-dir ... --benchmark bc4 bc5` | step 10 |

### Code versions

* `src/grendel-sensitivity` (8802f50) is used as it is, no patch. `git diff c621279 8802f50 -- grendel/` touches
  only four documentation/metadata files (geometry README/SOURCE.json/manifest.json, eventcalc PROVENANCE.md);
  the production code of BC5 is also unchanged since 461b7e5, and the template generator
  (`reco/exhad_generation.py`, `models/scalar/templates_exhad.py`) since 536e586. The scan's `code_digest` is
  therefore different (ae2b430e... here, 16eaa1e2... in the published run_metadata), the numbers are not.
* `particle`: the BC4 pool needs 0.26.2 (b-hadron lifetimes enter the weights). `00_setup.sh` installs it
  `--no-deps --target work/scalar/overlay/particle-0.26.2`; only step 01 puts it on `PYTHONPATH`. Everything else
  uses the venv's pinned 1.0.1 (the BC5 production of Oct 2026 and the scans used 1.0.1). The control scans prove
  the choice (below): with 0.26.2 the BC4 rows are the published ones to 4e-16.
* exHad: the pinned e41b07b instead of 7bb2308. Upstream history (`git log 7bb2308..e41b07b`, two commits): the
  `data/scalar` tables are unchanged, the code changes are the ALP threshold 1.911 -> 1.910 and how the string
  kernel installs its fragmentation pointer; the upstream golden digests of every scalar scenario (including
  `scalar-1809-had-2.30`, string fragmentation) are the same at both commits, only ALP digests changed. Our
  templates confirm it (step 06).
* exHad's venv has `particle` 1.1.0; exHad takes only identity data (charges) from it, and the bundle field that
  uses grendel's model through it (`ctau_m_u2eq1_reference`, the Winkler lifetime) is a diagnostic the v2.5 scans
  do not read.
* Pythia: the shower is `higgs_shower.cc`, compiled with the conda `c++` against `third_party/pythia8317`; the
  original used the same C++ program (no `pythia8mc` needed).

## Steps

| step | what | output (`work/scalar/`) | runtime here |
|---|---|---|---|
| `00_setup.sh` | folders, `particle` 0.26.2 overlay | `overlay/` | seconds |
| `01_bc4_vectors.sh` | BC4 four-vectors (`py/produce_bc4.py` = `scan.produce_missing`, atomic writes) | `bc4/llp_4vectors/` 86 x 122 MB | 15 min, 1 slot, 1 GB |
| `02_bc5_higgs_mg5.sh` | MG5 LO ggF/VBF/WH/ZH/ttH, 200 000 events each, `--nb-core 4` (as published) | `bc5/madgraph/` | 28 min, 4 slots |
| `03_bc5_higgs_shower.sh` | Pythia 8.317 shower, seed 1, of our LHE and (check) of the stored LHE | `bc5/madgraph/higgs_*.csv`, `bc5/shower_of_published_lhe/` | 2 x 12 min, 1 slot |
| `04_bc5_higgs_hist.sh` | `higgs_pool` to work/ (not into src/) | `bc5/higgs_hist/` | seconds |
| `05_bc5_vectors.sh` | BC5 four-vectors (`py/produce_bc5.py` = `scan_bc5.produce_missing`), tracked histogram | `bc5/llp_4vectors/` 112 x 181 MB | 30 min, 1 slot, 1.5 GB |
| `06_templates.sh` | exHad scalar-1809, 112 masses | `templates_exhad/` | 20 min on 4 slots |
| `07_scan_controls.sh` | control masses with the published and with our templates | `<bench>/scans/{control_reftpl,main}/` | see below |
| `08_scan_all.sh` | all 198 one-mass jobs (`scan_one.sh` -> `scan_job.sh`) | `<bench>/scans/main/out/<m>/` | see below |
| `09_collect.sh` | tables, nsig10, BC5 run_metadata + paper_numbers, comparisons (`py/compare_rows.py` per row, `py/compare_tables.py` per column) | `results/{bc4,bc5}/` | seconds |
| `10_figures.sh` | Figs. 4-5 from our and from the published tables, pixel (PNG) and byte (PDF) comparison; copies ours to `figures/scalar/` | `figures/` | seconds |
| `11_higgs_hist_impact.sh` | BC5 1, 15, 50 GeV with the tracked vs the rebuilt Higgs histogram (validation of the choice in 05) | `bc5/hist_impact/` | ~30 min, 1 slot per job |

`scan_one.sh <bc4|bc5> <mass> [RUN]` is one job (a per-mass lock prevents duplicates; finished masses are skipped);
`scan_job.sh` runs it inside its slot behind a memory gate (jobs start one at a time when >= 3 GB are free).

## Validation

All reports are in `work/scalar/validation/`. Tolerances are stated per check.

**Four-vectors.**
* BC5 (step 05, `bc5_vectors.txt`): the stored `asbuilt_20261014/bc5/llp_4vectors` are readable here. 4/112 files
  are byte-identical; in the other 108, 12 lines (median; max 28) of ~1.7 million differ, each by one unit in the
  9th significant digit (`%.8e`): largest difference 9.5e-9 (weight relative to itself, E/p relative to E).
  Same line count and row structure everywhere. This is libm last-bit noise (cosh/sinh/cos of macOS arm64 vs glibc
  x86_64) reaching the printed digits; tolerance 1e-7 and < 1e-4 of the lines: PASS.
* BC4 (step 01, `bc4_vectors.txt`): the stored pool files are HFS-compressed on the Mac drive and read as empty
  here, so only their sizes are known (`v2_3_20261016/lxplus/inputs_manifest.txt`): 86/86 sizes equal (row
  counts and the sign of every number, i.e. the RNG stream and the kinematics). The contents (weights) are
  checked by the control scans: with `particle` 0.26.2 the rows are the published ones to 4e-16 (below); with
  1.0.1 the weights would be 0.1-0.4 % lower.
* Per-row input digests of the published `run_metadata.json` (SHA-256 of vectors + template) cannot match: both
  files differ in the last bits. The digest method itself was checked: the stored BC5 vectors + stored template
  reproduce the published digest of BC5 0.14 GeV.

**Higgs histogram (steps 02-04).**
* MG5 (`mg5.txt`): every fresh process folder gets MG5's seed offset 21, as published. The LO cross sections are
  the published ones to all printed digits for ggF (18.9424 pb), WH, ZH and ttH; VBF gives 4.72139 pb against
  4.72422 pb (4.721 +- 0.003 vs 4.724 +- 0.004 pb: the same within the integration error). Event by event: ggF,
  WH, ZH and ttH are the published events (60-75 % byte-identical, the rest with momenta equal to <= 7e-5 of the
  particle energy, the same particles); VBF (2 -> 3, many channels) adapts a different grid and gives a
  statistically equivalent, different sample.
* Shower (`shower_of_published_lhe.txt`): the stored 2026-09-25 LHE files showered here (Linux) against the
  published shower of the same files (macOS). WH is the published one (36 of 200 000 rows differ, none by more
  than 1e-6 relative). In ggF, VBF, ZH and ttH 73-99 % of the rows differ: once a last-bit difference flips one
  accept/reject decision of the shower, the history (and the random-number sequence) of that event diverges, so
  these are other shower realisations of the same hard events. Mean Higgs pT: 51.729 / 51.714 GeV (ggF),
  93.857 / 93.819 (VBF), 89.551 / 89.535 (ZH), 130.164 / 130.082 (ttH), i.e. within 0.02-0.06 %.
* Histogram (`higgs_hist.txt`, two-sample chi^2 with cells merged to >= 10 counts): (b) those showers binned as the
  tracked `higgs_pt_y_14TeV.csv.gz`: WH identical cell by cell; ggF, VBF, ZH, ttH chi^2/ndf 0.63-0.84 (below 1
  because both samples share their hard events), binned mean pT within 0.09 %, mean |y| within 0.2 %. (a) the whole
  chain rebuilt here (02 MG5 + 03 shower + 04): chi^2/ndf 0.61-0.94 per mode (p >= 0.965); binned mean pT within
  0.35 % (ggF 51.70 vs 51.88 GeV, VBF 94.38 vs 94.17 GeV), mean |y| within 0.51 % (VBF 1.183 vs 1.177).
  Statistically the same histogram, not the same numbers.
* Which histogram the BC5 four-vectors use: the tracked one (05's default `HIGGS_HIST=tracked`), as the published
  production did. It is the published input itself, which this machine cannot regenerate bit for bit (shower
  divergence above), and only with it can the BC5 rows be checked against the published ones to the last bit
  (step 07: 15 and 50 GeV bit-identical). Step 11 (`higgs_hist_impact.txt`) measures what the rebuilt histogram
  would change, with the same RNG stream (seed 42, masses 1, 15, 50 in their own stream) and our templates:
  1 GeV u2_min x1.0028, u2_max x1.0034, peak N x1.0059; 15 GeV x0.9909, x1.0202, x1.0264; 50 GeV x1.0115,
  x0.9970, x0.9730; N >= 10 edges within 1.7 %; every row sensitive at N >= 3 and N >= 10 with either. That is the
  scan's own Monte Carlo precision: these rows rest on 2 900-4 500 h -> SS hits (effective sample size at the
  peak 1 900-3 200, 1.8-2.3 %), which the two histograms place differently.

**exHad templates (step 06, `templates.tsv`).** All 112 bundles: mass, exHad lifetime `ctau_m_u2eq1`, n_templates,
base seed and the per-point exHad seed equal to the published bundles. No bundle is byte-identical, but in 106 of
112 every one of the 20 000 events has the same daughters with momenta equal to 1e-6 (the same decays up to
floating-point noise). The other 6 differ in a few events or one sub-request: 2.4 and 3.1 GeV (the published
bundles retried a sub-request after exHad's four-momentum check), 2.5 and 3.6 GeV (ours retried one), 4.2 and
11 GeV. Statistically: channel fractions, mean multiplicity, mean charged multiplicity and P(>= 2 charged)
agree with |z| <= 0.56 everywhere; two-sample KS of the daughter |p| and charged-daughter |p|: p >= 0.97.

**Control scans (step 07, `controls_*.txt`).** BC4 0.5, 1, 2, 3.6 GeV and BC5 0.3, 1, 15, 50 GeV at full
statistics:
* with the published templates (`control_reftpl`): BC4 every edge column (u2_min, u2_max, N >= 10 edges)
  bit-identical, all other columns within 1.6e-15 (peak N, effective-sample-size diagnostics: summation order);
  BC5 15 and 50 GeV rows bit-identical, 0.3 and 1 GeV within 5e-15.
* with our templates (`main`, the production rows): identical to the above except BC4 3.6 GeV (our bundle retried
  one sub-request): u2_min x0.9977, u2_max x1.0031, peak N x1.0019 (3.73 -> 3.74; still sensitive, as published).

**Final tables (step 09).** All 198 one-mass jobs of 08 finished at the first attempt (no failure, no rerun). Reports:
`final_{bc4,bc5}.txt` (per row, `py/compare_rows.py`), `final_{bc4,bc5}_columns.txt` (per column,
`py/compare_tables.py`), `final_{bc4,bc5}_nsig10.txt`, `final_bc5_run_metadata.txt`, `sec4_numbers.txt`.
* No verdict changes anywhere. BC4: 75/86 masses sensitive at N >= 3 and 69 at N >= 10, the same masses as
  published; BC5: 112/112 at both, and all 14 `*has_sensitivity*` columns (channels, BR(h->SS) = 1e-3 overlay)
  equal row by row.
* The rows split exactly along the decay templates (step 06). At the 81 BC4 and 106 BC5 masses whose bundles are
  event by event the published ones, every edge column is bit-identical (BC4 10, BC5 70 columns: u2_min, u2_max,
  peak_u2, their N >= 10, open-edge and per-channel/overlay versions), as are n_events, n_hits, ctau and the flags;
  only peak N, the per-channel yields at the edges and the effective-sample-size diagnostics differ, by at most
  8.7e-15 (BC4) and 1.6e-14 (BC5) relative (floating-point summation order). 16/86 BC4 and 43/112 BC5 rows are
  bit-identical in every column; no table is bit-identical as a file.
* The masses whose bundles differ (2.4, 2.5, 3.1, 3.6, 4.2 GeV, and 11 GeV in BC5) move by the decay Monte
  Carlo: BC4 u2_min x0.9977-1.0011, u2_max x0.9993-1.0031, peak N x0.9968-1.0019 (N >= 10 edges within 0.18 %);
  BC5 u2_min x0.9970-1.0021, u2_max x1.0005-1.0110 (the largest at 4.2 GeV), peak N x0.9982-1.0035. Diagnostics
  of these rows (sub-channel yields at an edge, largest single-event fractions) move by up to 41 %, as expected
  for quantities carried by a handful of events.
* `sensitivity_nsig10.csv` (the N >= 10 projection by `grendel.io.thresholds`): against the published files the
  differences reach 1.9e-13 (BC4) and 7.2e-13 (BC5) relative also in the template-equal rows. That is the parser,
  not the data: `thresholds` re-reads the table with pandas' default float parser, which on the published machine
  (pandas 3.0.6, macOS arm64) drops the last digits (`0.0011888522982516816` -> `0.0011888522982516`), while the
  pinned pandas 3.1.1 here parses exactly. Projected here, the published `sensitivity.csv` gives a nsig10 table
  that ours matches exactly like the main tables (8.7e-15 / 1.6e-14 outside the template masses). pandas 3.0.6
  installed on this machine reproduces 1759/1827 (BC4) and 8832/9153 (BC5) published nsig10 cells, the rest
  differing in the 17th significant digit (presumably x86_64 vs arm64 arithmetic in the parser), so these two
  files are not reproducible bit for bit here.
* Sec. IV numbers. BC5 `paper_numbers.txt` equals the published one except one line: MATHUSLA40 "GRENDEL upper edge
  4.8-13.8x larger" (published 13.7x): the maximum sits at 11 GeV, whose u2_max is x1.0029 (factor 13.728 ->
  13.768). The paper text does not quote that number. BC4 `edge_ratios.py`: identical except the 3.0-4.0 GeV
  FASER2 lower-edge minimum 6.27 (6.25) and the SHiP upper-edge maximum 67.09 (67.30), both from 3.6 GeV. Every
  number quoted in Sec. IV is reproduced: BC4 deepest 7.2e-12 at 1.0 GeV, 1.4e-11 over 1-3.5 GeV, reach 3.6 GeV,
  FASER2 lower edge 3.1-9.7x higher, SHiP 0.08-0.52; BC5 2.05e-17 at 50 GeV, peak N >= 95 at every mass,
  SHiP below 2.5 GeV 2.1-8.9x smaller (the text says 2-8; the published `paper_numbers.txt` has the same
  2.1-8.9), h->SS bands over 2-55 GeV.
* BC5 `run_metadata.json`: 26 of 33 fields equal; the others are the code digest (ae2b430e... vs 16eaa1e2..., see
  Code versions), the three input paths (`work/scalar/...` vs `inputs/...`), `workers` (named after this
  pipeline), `total_time_s` (37 888 s here, 84 305 s on lxplus) and the 112 per-mass input digests (none equal:
  vectors and templates differ in the last bits, see above).

**Figures (step 10).** The renderer here, fed the published v2.5 tables, gives PNGs pixel-identical to the
Overleaf `bc4_grendel_paper.png` and `bc5_grendel_paper.png` (0 differing pixels; the files differ only in PNG
metadata). Our tables, rendered the same way (`work/scalar/figures/ours`, copied to
`figures/scalar/`): BC4 0.059 % of the pixels differ from the Overleaf PNG (largest channel difference 19/255, none
above 32), all in x 1216-1463, y 470-650, i.e. on the contour lines between 2.4 and 4.2 GeV; BC5 0.068 % (largest
36/255, 0.0002 % above 32), on the contour lines at 2.4-4.2 and 11 GeV. Our PDFs differ from the Overleaf PDFs in
content (36 146 vs 36 144 and 34 927 vs 34 921 bytes without the date stamp), while the PDFs of the published
tables rendered here are identical to them apart from `/CreationDate`.
Attribution: our tables with only the rows of the template-differing masses replaced by the published rows
(`py/splice_rows.py`, `work/scalar/figures/spliced`) render pixel-identical to both paper PNGs (0 differing pixels),
so every visible difference comes from the decay Monte Carlo of those five (BC4) / six (BC5) masses.

## Runtimes (this laptop: Ryzen 5 PRO 7530U, 6 cores / 12 threads, shared with four other pipelines)

| step | measured |
|---|---|
| 01 BC4 four-vectors | 15 min (1 slot; I/O bound, 13 s per 122 MB file) |
| 05 BC5 four-vectors | 30 min (1 slot) |
| 06 templates | 21 min on 4 slots (1-2 s per mass below 2 GeV, 20-35 s for string fragmentation above) |
| 02 MG5 | 28 min on 4 slots (ggF 1, VBF 19, WH 2, ZH 4, ttH 3 min; published run on macOS: about the same) |
| 03 shower | 12 min per pass of the five modes (1 slot; 1.6-3.5 min per mode), two passes |
| BC4 scan, one mass | 5.0-12.1 min (median 7.1), 1.4-2.3 GB RSS (lxplus: 13-27 min) |
| BC5 scan, one mass | below 5 GeV 1.9-11.2 min (median 6.8), 0.8-2.4 GB; 5-62 GeV 1.2-3.3 min, 0.5-0.8 GB (lxplus: 17-29 min) |
| 08 all 198 scans | 21.7 slot-hours (BC4 11.0, BC5 10.7; lxplus 71.4 core-hours); 14.2 h elapsed (2026-11-06 18:06 to 11-07 08:21) with 0-8 of the 11 shared slots (0-2 from about 20:30 to 00:00, while long FONLL, HNL and SensCalc jobs held the rest). Alone on the machine (estimate): SCALAR_JOBS=6 and the 3 GB memory gate allow 4-6 jobs at once, so about 4-5 h |
| 11 histogram impact | 6 scans, 1.5-6.2 min each, plus two 3-mass productions |


## Discrepancies and notes

* **The rebuilt tables are not the published files bit for bit, and cannot be on this machine.** Everything
  upstream of the decay templates reproduces the published numbers: at the 81 BC4 and 106 BC5 masses with
  equivalent templates every edge is bit-identical. The remaining differences come from exHad on Linux x86_64 vs
  macOS arm64: same seeds, but floating-point noise changes a few events or one retried sub-request in 6 of the
  112 bundles (2.4, 2.5, 3.1, 3.6, 4.2, 11 GeV), which moves the edges of those rows by up to 0.6 % (u2_max at
  4.2 GeV: 1.1 %) and their peak N by up to 0.35 %. No verdict, no number quoted in the paper and no figure pixel
  outside those masses changes (step 10 attribution).
* **BC4 four-vectors.** The stored BC4 pool is unreadable here (HFS compression lost on copy), so step 01 is
  checked by file sizes only; the final table closes that gap: the bit-identical edges at 81 of 86 masses need the
  same weights to the last bit, which only the `particle` 0.26.2 overlay gives (with the pinned 1.0.1 every BC4
  row would move by 0.1-0.4 %).
* **Higgs histogram.** The BC5 production uses the tracked histogram (the published input). The one rebuilt by
  02-04 is statistically equivalent (VBF: MG5 adapts another integration grid here, 4.72139 vs 4.72422 pb, a
  different sample; the other modes share 60-75 % of their events byte for byte with the published LHE; the
  Pythia shower diverges event by event between the two platforms). Using it would move the BC5 rows by up to
  2.7 % (step 11), within the scan's own Monte Carlo precision.
* **`sensitivity_nsig10.csv`** inherits a pandas float-parser difference between the published machine and this
  one (up to 7.2e-13 relative); see Final tables.
* **Paper text vs the published pipeline output** (not a reproduction issue, noted for the authors): Sec. IV says
  SHiP reaches mixings "2--8 times smaller" below 2.5 GeV; `paper_numbers.txt` (published and ours) gives 2.1-8.9.
* **How it was run.** The steps were started one by one with the scripts of this folder (07 with SCALAR_JOBS=3,
  08 with 5, raised to 8 while running), 11 alongside 08; `run_all.sh` runs the same scripts in the same
  dependency order. No result depends on the parallelism: every scan job, exHad point and sub-request has its own
  seed, and the four-vector productions run their whole RNG stream in one process.
