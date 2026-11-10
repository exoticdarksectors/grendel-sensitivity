# FONLL: the heavy-flavour production grids

The BC4, BC5, BC10 and HNL production samplers read two FONLL (pT, y) grids, `pT y dsigma/dpT/dy` in pb/GeV for
pp at 14 TeV with `NNPDF40_nlo_as_01180`: B hadrons (bottom) and the D0-convention charm grid. They are tracked
in the code (`src/grendel-sensitivity/grendel/production/data/fonll/central/`) and the other areas use those
files. This area proves that they can be regenerated from scratch on this machine, with FONLL 1.3.3 and the
patches the code ships (`grendel/production/fonll_grids`).

Run it all (resume-safe, CPU through `slot_run`, at most 4 slots):

    source env/activate.sh
    setsid nohup bash pipelines/fonll/run_all.sh > logs/fonll/run_all.log 2>&1 < /dev/null &
    source pipelines/fonll/env.sh && python pipelines/fonll/progress.py       # while it runs
    source pipelines/fonll/env.sh && python pipelines/fonll/05_compare.py     # any time: compares what exists

**Reference material.** The original FONLL workspace (`grendel/shared/NNPDF40/fonll-local`: raw quark grids,
fragmentation outputs, calibration) is read from the optional reference tree `$GRENDEL_REF`
(`pipelines/common/ref.sh`, `reference/README.md`; `FONLL_ORIG_WS` overrides it). Without it, `xcheck_fragment.py`
and sections 2-4 of `05_compare.py` print a SKIP line and the run still finishes; the comparison with the tracked
grids in `src/grendel-sensitivity` (section 1) and the CTEQ6.6 numbers always run. No step uses it as an input.

## What the published production reads

* `grendel.production.fonll.fonll_parser` resolves the grid of each quark as `GRENDEL_FONLL_BOTTOM_GRID` /
  `GRENDEL_FONLL_CHARM_GRID` if set, else the tracked
  `data/fonll/central/fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_{bottom,charm}.dat`.
  The variables are only used by the uncertainty-band campaigns (`grendel.band.*`); no central production
  recipe of the v2 tables sets them.
* The grids were added to the repo once (ccf7ea0, 2026-09-21) and never changed up to 8802f50. The generator
  changed afterwards only in non-numerical ways (resume logic, `--calibrate-only`, manifest notes, the Linux
  make flags of `install.py`); `PINS.json` changed only its upstream URL.
* Every stored copy in the archive has the same data rows: 80 files (40 per quark) in the reference tree —
  `fonll-local/output`, `fonll-nnpdf40/grids` (v0.1.0, eb5564b, 2026-06-01), the model repos the v1 and v2 pools
  were produced with (`signal-models/*/hnl/data/production/fonll/central`, `grendel-mono`, the
  `llpatcolliders_*` worktrees), and every campaign code snapshot (`v2_baseline_20260929/code`,
  `v2_3_20261016/code`, `hnl_histat_20261006/code_main`, `bc4_winkler_v2_20261007/code_main`, ...). SHA-256 of the
  non-`#` lines: bottom `ec33ab8bc3e5...`, charm `feb64ed41665...` (full values in `05_compare.py`). The files
  differ only in three header lines (scale wording, the `lhapdf_member`/`variation_*` lines of the later
  generator, and a hand-edited `meson_mass` note in the tracked copy).
* They were computed on 2026-05-31 in `grendel/shared/NNPDF40/fonll-local` (macOS arm64, Homebrew
  gfortran 15.2.0 at `-O`, LHAPDF 6.5.6) with `scripts/generate_meson_grids.py --pdf nlo --quark {bottom,charm}
  --grid-workers 4`, the script `grendel.production.fonll_grids.generate` is a port of (identical input cards,
  grids, fragmentation cards, feeddown combination and fit; checked function by function). Its raw quark grids
  and fragmentation outputs are stored in `fonll-local/run/nlo_{bottom,charm}/` and are the intermediate
  references below. A 1-process rerun there on 2026-06-10 gave byte-identical raw grids, so the rapidity
  chunking does not change the numbers.

## Steps

| Step | Does | Inputs | Outputs (`work/fonll/`) | Runtime here |
|---|---|---|---|---|
| `env.sh` | paths sourced by every step: `GRENDEL_FONLL_DIR`, `GRENDEL_FONLL_GRID_DIR`, `LHAPDF_DATA_PATH=work/lhapdf:work/fonll/lhapdf` | | | |
| `00_fetch_pdfs.sh` | `cteq66` from the LHAPDF server; SHA-256 of `cteq66{.info,_0000.dat}` and `NNPDF40_nlo_as_01180{.info,_0000.dat}` checked against the files the original workspace read | network | `lhapdf/cteq66` | seconds |
| `01_build_fonll.sh` | copies `third_party/fonll` (c7086e4, kept pristine), applies the 4 shipped patches with `install.py` (hash-checked before/after) and builds `fonllgridlha`, `fragmfonll` at `-O` (see Build notes) | conda gfortran 15.3, LHAPDF 6.5.6 | `fonll-1.3.3/Linux/` | 25 s |
| `smoke.py` | the whole chain (raw grid, calibration + fit, both fragmentation passes, D* feeddown, writer) on a 4 y x 5 pT grid with 20 raw pT nodes | 01 | `smoke/` | 3 min (bottom), 8 min (charm) |
| `02_calibrate.sh` | unit `calib`: charm CTEQ6.6 D* feeddown calibration (`generate --calibrate-only`); unit `closure`: bottom CTEQ6.6 public point (`closure.py`) | 00, 01 | `grids/output/charm_feeddown_calibration.json`, `validation/bottom_cteq66_public_point.json` | calib 20 min, closure 17 min (machine full) |
| `03_raw_quark_grids.sh` | NNPDF4.0 raw quark grids, 100 y x 80 pT per quark, as 4 rapidity chunks per quark, each its own 1-slot process (`raw_chunk.py`), then merged | 01 | `grids/run/nlo_{bottom,charm}_central/nl_{b,c}_.out` | 3.5 min per rapidity row at central y, 6.1 min near abs(y) = 3 (machine full): 1.5-2.5 h per chunk, ~8 single-core h per quark |
| `04_fragment_central.sh` | `generate --pdf nlo --quark {bottom,charm} --reuse-existing-grids`: fragmentation and grid writer on the step-03 raw grids | 02 (charm), 03 | `grids/output/fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_{bottom,charm}.dat` | bottom 30 s, charm 2 min |
| `05_compare.py` | all the comparisons below (runs on whatever exists) | | `validation/compare.json`, `validation/summary.txt` | seconds |
| `xcheck_fragment.py` | validation: the stored macOS raw grids (`fonll-local/run/nlo_*/nl_*_.out`) and calibration fragmented here with the shipped `generate --reuse-existing-grids` | 01 | `xcheck/` | 3 min |
| `progress.py` | while step 03 runs: points done, ETA, and the running comparison of the rows done with the macOS raw grid | | stdout | seconds |
| `run_all.sh` | 00, 01, one queue of 1-slot units (calib, closure, 8 chunks) at most 4 at a time, merges, 04, `xcheck_fragment.py` (once), 05 | | `logs/fonll/run_all.log` | ~17 slot-hours; ~4.5-5 h wall with 4 slots |

`task.sh <unit>` runs one unit (`calib`, `closure`, `raw:<quark>:<i>`, `merge:<quark>`, `frag:<quark>`) under one
slot and a per-unit flock; every unit skips itself when its output is complete, so any step or `run_all.sh` can
be restarted at any time. The upstream driver (`generate --grid-workers N`) runs the N rapidity chunks as threads
of one process holding N slots and restarts all of them when interrupted; `raw_chunk.py` makes the same chunks
(same cards, same merge) separate 1-slot processes, which the shared slot pool schedules better and which resume
individually. The rest is the shipped generator unchanged.

Logs: `logs/fonll/<step>.log`, `logs/fonll/run_all.log`; per-FONLL-command logs in `work/fonll/grids/logs/`.

## Build notes

* The installer's Linux path passes `F77=gfortran -Wall -fno-automatic -fallow-argument-mismatch SYSOBJ=`, but
  the Linux block of the FONLL Makefile leaves `OPT` empty, i.e. `-O0 -ggdb`, while the original macOS build
  used `-O` (Darwin block). `01_build_fonll.sh` adds `OPT=-O` to the installer's make command (a 3-line wrapper
  around `install.py`; nothing in `src/` is modified; `FONLL_OPT=` gives the bare installer). Test (`buildtest_rows.py`, `buildtest_cmp.py`) on two
  rapidity rows (160 raw-grid points) of the bottom grid against the stored macOS raw grid:
  * `-O`: matched FONLL cross section within 1.6e-5 (median 1e-11) of the macOS value; the massive and
    massless fixed-order columns agree to 3e-9 and 1.7e-6, the resummed column to 1.2e-5;
  * `-O -mfma -ffp-contract=fast` (to mimic arm64's fused multiply-add): byte-identical to `-O` on the 100 points compared;
  * `-O0` (bare installer): ~7x slower, and the massless column already off by up to 1.4e-4 on the 16 points it
    computed before it was stopped.
* Raw grids are not byte-identical to macOS: arm64 vs x86-64 and Apple libm vs glibc differ in last bits, and
  FONLL's adaptive integration propagates them. The differences are far below FONLL's own integration error
  (quoted in the raw grid as `ehrs`).
* Speed: one raw-grid rapidity row (80 pT points) takes 3.5 min (central y) to 6.1 min (abs(y) near 3) on one
  slot here with the machine full (the M-series Mac: ~1.6 min), so a grid is ~8 single-core hours, not the 2 h of
  the upstream README.

## Validation

`05_compare.py` writes `work/fonll/validation/summary.txt` and `compare.json`. Final results of the full run
(2026-11-06 17:41 -> 2026-11-07 01:32 on the 12-core laptop, sharing the CPU with the other pipelines):

**Result: the regenerated central grids reproduce the tracked ones to 3.4e-7 (bottom) and 3.1e-5 (charm) per bin,
and the integrated cross section to 6e-11 (bottom) and 5e-8 (charm).** The data rows are not byte-identical; the
only source of the differences is the last-bit arithmetic of FONLL's numerical integration on x86-64/glibc vs
arm64/Apple (the fragmentation and the writer are exact, see the cross-check row).

| Check | Reference | Result |
|---|---|---|
| Inputs: `cteq66` and `NNPDF40_nlo_as_01180` (`.info`, member 0) | files read by the original workspace | SHA-256 identical |
| Patched sources | `patches/PINS.json` | the 4 files have the pinned patched SHA-256 (as `fonll-local/src/fonll`) |
| Smoke (`smoke.py`) | | both quarks: full chain runs, grids written with the expected rows |
| Fragmentation half of the chain (`xcheck_fragment.py`): stored macOS raw grids + original weight, fragmented and written here | tracked grids | **data rows byte-identical**, 10000/10000 for bottom and for charm (sigma and every bin identical). `fragmfonll` outputs: dsigma column identical at all 9900 / 9900 / 29700 points; 1-3 % of the lines differ only in the last digit of the error column |
| Bottom public point (`closure`): CTEQ6.6, B hadron, pT = 5 GeV, y = 0 | original macOS run 7.5352771279e6 pb/GeV; public FONLL web 7.5329e6 | 7.5352771274e6: -6.6e-11 vs macOS; +0.032 % vs public (as recorded in `PROVENANCE.md`) |
| Bottom CTEQ6.6 raw rows y = -0.5, 0, 0.5 | `fonll-local/run/cteq66_bottom/ct_b_.out` | 119/240 lines byte-identical; matched cross section max rel. diff 2.0e-9 |
| Charm calibration raw grid (CTEQ6.6, 3 y) | `fonll-local/run/validate_cteq66_charm_y0/val_c_.out` | 105/240 lines byte-identical; matched cross section max rel. diff 3.5e-7 (1.7e-3 of FONLL's integration error) |
| Charm calibration fragmentation (9 direct + 27 D* points) | stored `fragmfonll_calibration_*.dat` | max rel. diff 2.1e-10 |
| Charm feeddown weight | 1.2152533641182024 (tracked header: 1.21525336412) | 1.2152533688784 (+3.9e-9); max residual vs the 9 public D0 points 1.9410e-3 (same) |
| NNPDF4.0 bottom raw grid, 8000 points | `fonll-local/run/nlo_bottom/nl_b_.out` | 2397/8000 lines byte-identical; matched cross section max rel. diff 6.1e-7, median 9e-12 (at most 5.0 x FONLL's own integration error at one point) |
| NNPDF4.0 charm raw grid, 8000 points | `fonll-local/run/nlo_charm/nl_c_.out` | 2068/8000 lines byte-identical; matched cross section max rel. diff 1.2e-4, median 7e-11 (at most 8.6 x the integration error). Single columns differ more where they nearly cancel (hdml 5.2e-2, ehrs 6.7e-2 max), not in the matched sum |
| `fragmfonll` outputs on those raw grids | stored macOS `fragmfonll_*.dat` | bottom central max rel. diff 3.4e-7; charm direct 2.3e-5, D* vector 3.5e-5 |
| **Final bottom grid** (B hadrons) | tracked `..._central_bottom.dat` | 4323/10000 data rows identical; per bin max 3.4e-7 at (pT, y) = (9.6, -2.88), p99 9e-9, median 2e-11; sigma(pT 0-50, abs(y) < 3) = 3.059034541e8 pb vs 3.059034542e8 (-6.1e-11); pT and abs(y) slices within 3.3e-10 |
| **Final charm grid** (D0 convention) | tracked `..._central_charm.dat` | 101/10000 data rows identical; per bin max 3.1e-5 at (pT, y) = (2.0, -2.76), p99 3.4e-7, median 1.5e-10; sigma = 2.722259249e9 pb vs 2.722259386e9 (-5.0e-8); slices within 5.9e-8 |
| Two-row `-O` build test (y = -3 and 0.03, bottom) | same | max 1.6e-5 on the matched cross section (one point); median 1e-11 |

"Matched cross section" is FONLL's `hdmv + G (hdrs - hdml)`, `G = pT^2 / (pT^2 + 25 m^2)`, the quantity
`fragmfonll` convolves (`xsecmtch`); the integration error is `G * ehrs` from the same row.

Why the final grids differ at all: since the fragmentation and the writer reproduce the tracked grids bit for
bit from the macOS raw grids, the regenerated grids can differ from the tracked ones only through the last-bit
differences of `fonllgridlha` (x86-64/glibc vs arm64/Apple libm), measured so far at <= 1.6e-5 per raw point
(median ~1e-11), and, for charm, through the 3.9e-9 shift of the fitted weight. Full-file hashes cannot match in
any case: the shipped generator writes three header lines differently from the tracked files (`scale:` wording,
the `lhapdf_member` / `variation_kind` / `variation_tag` lines, and the tracked file's hand-edited `meson_mass:`
note), and the charm header records the weight to 12 digits. The comparison is therefore on the data rows
(SHA-256 of the non-`#` lines), per bin, and on the integrated cross section the samplers use
(`get_sigma_total`), in total and in pT / abs(y) slices.

## Not reproduced here

* The scale / PDF-replica / mass / alpha_s variation campaign (`generate --campaign`, `combine`): the figures use
  the central grids only.
* The public FONLL web-form values (CTEQ6.6) are taken as recorded in the code (queried 2026-05-29); they are not
  re-queried.

## Run record and resume

`run_all.sh` ran on 2026-11-06 17:41 -> 2026-11-07 01:32 (`logs/fonll/run_all.log`): 7 h 51 min wall with at
most 4 slots; the 8 raw-grid chunks took 2.0-3.9 h each (bottom 7288-9998 s, charm 11170-14039 s with the machine
full). Check a run with

    source env/activate.sh && source pipelines/fonll/env.sh && python pipelines/fonll/progress.py
    grep -E 'start|finished' logs/fonll/run_all.log
    cat work/fonll/validation/summary.txt          # after step 05

Resume after any interruption by launching `run_all.sh` again: finished units are skipped, an interrupted
rapidity chunk restarts from its beginning. Do not edit `task.sh` or `run_all.sh` in place while they run
(bash reads scripts as it executes them).

Using the regenerated grids downstream (not needed: the other areas use the tracked files, which these
reproduce): `GRENDEL_FONLL_BOTTOM_GRID` / `GRENDEL_FONLL_CHARM_GRID` = the two files in `work/fonll/grids/output/`.
