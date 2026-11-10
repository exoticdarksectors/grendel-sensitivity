# pipelines/senscalc — Wolfram Engine and SensCalc

Everything behind the paper figures that was computed with Mathematica, rerun with the free Wolfram Engine on Linux.

| Mathematica product | Used in | Published recipe (reference, read-only) | Step here |
|---|---|---|---|
| FASER2-FPF BC10 curve `curves/data/bc10/projections/faser2_senscalc_fpf.dat` (SensCalc eb2c456, headless) | Fig. 7 | `grendel/campaigns/bc10/faser2_senscalc_20261025` (`README.md`, `scripts/`, `logs/`, `out/`) | 11–14 |
| SHiP BC10 curve `curves/data/bc10/projections/ship_senscalc.dat` (M. Ovchynnikov's JSON, 2026-10-16; private communication) | Fig. 7 | the same campaign reran it with the same pipeline: both edges within 0.5 % (17 masses) | 11–14 |
| `grendel/models/alp/data/senscalc_2501/*` (decoded SensCalc v1.4 MX decay tables; Γ(a→μμ) reference of the BC10 exHad templates) | Fig. 7 | `PROVENANCE.md` there; `export_senscalc_2501.py` + `.wls`, commit db79c0e (unchanged up to 8802f50) | 20 |
| SensCalc v1.4 LHC production inputs (`export_senscalc_2501_production.py`; BC10 production audit, nothing tracked) | – | same tool | 21 |
| `grendel/models/scalar/data/senscalc_quartic/*` (B→X_s SS, B_s→SS, B⁺→K⁺SS) | Fig. 5 | `PROVENANCE.md` there: upstream data files copied verbatim | 22 |
| (d) `grendel/models/scalar/data/eventcalc_1809/ctau.json.gz` (exHad scalar-1809 = EventCalc/SensCalc lifetime) | Figs. 4, 5 | `PROVENANCE.md` there (copied from exHad); origin SensCalc `ctauS.mx` | 23 (cross-check) |

**Reference material.** The campaign (`grendel/campaigns/bc10/faser2_senscalc_20261025`, its SensCalc tree with every
output and its logs), the stored July 2026 production export and M. Ovchynnikov's SHiP JSON are read from the optional
reference tree `$GRENDEL_REF` (`pipelines/common/ref.sh`, `reference/README.md`). Without it every comparison with them
prints a SKIP line (in `work/senscalc/validation/summary.txt` and `work/skipped_validations.txt`) and the steps still
finish; the comparisons with the files tracked in `src/grendel-sensitivity` always run. The SHiP JSON is a private
communication, not public: it is the one input of a figure curve that only the reference tree holds. With it, step 14
regenerates `ship_senscalc.dat` (byte-identical to the tracked file); without it, step 14 copies the tracked file.

Other grep hits for SensCalc/Mathematica are not inputs: `curves/tools/vector_bc10.py` only says that no NA62-dump BC10 curve
is drawn because SensCalc needs Mathematica; `vector_codexb_hnl.py` reads a PDF made by Mathematica;
`curves.tools.senscalc_bc10_ship` converts the SHiP JSON (Python, step 14).

## Setup (no kernel needed)

| Script | What | Runtime | Status |
|---|---|---|---|
| `00_install_wolfram.sh` | Downloads Wolfram Engine **15.0.0** for Linux (the version of the reference runs; 2.3 GB, `WolframEngine_15_LIN.sh`, sha256 `5352d1e5…32ee`) and installs it in user space under `env/wolfram/WolframEngine` (8.2 GB). The installer refuses non-root installs unless given `-auto -wolframscript=n`; WolframScript 1.14.0 is unpacked from the bundled `.deb` with `dpkg-deb -x` into `env/wolfram/WolframScript`, and `env/wolfram/bin/wolframscript` is a wrapper that sets `WOLFRAMSCRIPT_KERNELPATH` to our kernel. `env/activate.sh` already puts `env/wolfram/bin` on PATH (also `WolframKernel`). | 2 min download + 3 min install | done |
| `01_fetch_sources.sh` | SensCalc (github.com/maksymovchynnikov/SensCalc) shallow fetches: `third_party/SensCalc` @ **eb2c456** and the worktree `third_party/SensCalc-v1.4` @ **4f7210b** (tag v1.4; `export_senscalc_2501.py` checks `git rev-parse HEAD`). FeynCalc **10.2.1** (tag `Release-10_2_1`, tarball sha256 `3594d33b…996f`) into `env/wolfram/Applications/FeynCalc` (SensCalc notebooks 1 and 2 do ``Get["FeynCalc`"]``). | 3 min | done |
| *activation* | `source env/activate.sh && wolframscript -activate` once, with the user's Wolfram ID. Writes `~/.WolframEngine/Licensing/mathpass` (outside the repo). | – | done 2026-11-06 |
| `02_check_engine.sh` | Version, licence limits (`$MaxLicenseProcesses`, `$MaxLicenseSubprocesses`), C compilation with `/usr/bin/gcc`, FeynCalc load. | 22 s | done: 15.0.0 Linux-x86-64, 2 kernels / 8 subkernels, gcc compiles, FeynCalc 10.2.1 loads in 14 s |

Validation of the setup:
- The eb2c456 checkout equals the GitHub tarball the campaign extracted (`src/SensCalc-eb2c456….tar.gz`, sha256 `de29db62…d705`):
  865 files compared by git blob hash, 0 mismatches (`py/compare_tarball_tree.py`).
- The v1.4 inputs match the SHA-256 the repository pins: the decay tables (widths, BRs, matrix elements, acceptance
  notebook) in `export_senscalc_2501.py`, the production inputs (coefficients, fragmentation, light-meson decays,
  bremsstrahlung, Drell-Yan, SM spectra, notebooks) in `export_senscalc_2501_production.py`; both `--check-only` pass.

## SensCalc headless (steps 10–14)

The campaign's method, unchanged: `wl/senscalc_headless.wl` evaluates the master notebooks cell by cell in a plain kernel,
emulates the front-end calls, and answers every dialog from preset tables in the drivers `wl/run_nb{1,2,3}_*.wls`
(see the campaign README, "Method"). The `wl/` files are copies of the campaign `scripts/`
(`wl/CAMPAIGN_ORIGINALS.sha256` holds the hashes of the copies as taken). Changes, all environment-only:
- `run_nb{1,2,3}_*.wls`, `check_soft.wls`: load `wl/senscalc_env.wl` before `LaunchKernels` and call
  `` SCEnv`subkernelEnv[] `` after it. That file puts `env/wolfram/Applications` (FeynCalc) on `$Path` of the main kernel and
  the subkernels and selects `/usr/bin/gcc` for `CompilationTarget -> "C"` (on macOS: FeynCalc in `$UserBaseDirectory`, clang).
- **Case-sensitive file system** (`py/case_links.py`, applied by `make_tree` in `lib.sh`). Two SensCalc notebooks build a
  path with `"production probabilities"` while the directory is `phenomenology/<model>/Production probabilities`; one is
  `codes/LLP distribution/prod-pheno-bremsstrahlung-proton.nb`, which defines the ALP-nucleon coupling
  `gppA[mLLP_, model] = Import[…/phenomenology/<model>/production probabilities/Coefficients-<model>-Lambda=1.-TeV.mx][[-1, 2]]`.
  APFS (macOS) is case-insensitive, ext4 is not: on Linux the Import returns `$Failed` without a message in the log,
  `Abs[gppA]^2` drops out of `FormFactorSquaredBrem`, and notebook 2 writes a bremsstrahlung production probability
  (`ProdProb_ALP-fermion_Bremsstrahlung-AP_<facility>.m`) too small by |gppA|² (×29 at 0.2 GeV, ×704 at 1 GeV; the
  missing factor is 10.7527² × the π⁰/η/η′-mixing polynomial with the eb2c456 couplings). Effect in the first Linux run
  (2026-11-06): bremsstrahlung N_ev ×0.001–0.03, FASER2-FPF upper edge −17 % at 1 GeV, SHiP-ECN3 upper edge −13 %
  at 1.2 GeV, lower edges ~1 % (B decays dominate there); outputs kept in `work/senscalc/diag_case_bug_20261107/`. Fix:
  a symlink `production probabilities -> Production probabilities` in every `phenomenology/<model>/` of the working tree
  (no SensCalc file changed; nothing lists those directories). The scan finds no other case-mismatched path read by the
  BC10 chain. The one-letter fix of both notebooks was merged upstream on 2026-11-07
  (maksymovchynnikov/SensCalc#4, merge commit e9a5a36). This pipeline keeps eb2c456, the version behind the
  published curves, and with it the symlinks; a SensCalc pin at or after e9a5a36 would not need them.
- `py/make_curve_files.py`: its glob `…Nev=3._Npot=*.json` also matches `*.metadata.json`; the campaign got the right file
  from the macOS listing order, ext4 lists the other one first. The copy skips `*.metadata.json`. With that fix it rewrites
  both campaign `.dat` files byte-identically from the campaign JSON, including the tracked `faser2_senscalc_fpf.dat`
  (sha256 `fa163868…f753`).

All runs: `$SENSCALC_SUBKERNELS` (default 4, as the campaign) subkernels, one run at a time through
`slot_run.py -n 4`, in a pristine `git archive` copy of eb2c456 under `work/senscalc/run/SensCalc@eb2c456` (SensCalc
writes its outputs into its own tree). Resume: markers in `work/senscalc/markers/`; delete one to redo that run.

| Script | Runs; runtime here (4 subkernels, `$ProcessorCount` = 1) vs campaign on macOS M-series | Output (under `work/senscalc/run/SensCalc@eb2c456`) |
|---|---|---|
| `10_smoke.sh` | separate tree `work/senscalc/smoke`: nb1 FASER2-FPF, nb2 LHC at 2·10⁴ events/channel (446 s, dominated by the deterministic bremsstrahlung integration), nb3 FASER2-FPF (891 s); criterion `py/smoke_check.py`: the chain runs and lands in the right place (deterministic channels exact, curve within 50 %, reach within 5 %; see Results) | – |
| `11_nb1_acceptances.sh` | nb1 SHiP-ECN3 cuts `ship` (1116 s; Mac 443 s), FASER2-FPF `none` (489 s; 169 s), FASER2-FPF-Etot100 `etot100` (763 s; 287 s); `inspect_acceptance.wls` metadata | `Acceptances/<exp>/…Default.m` |
| `12_nb2_distributions.sh` | nb2 LHC and SPS, B, Eta, EtaPr, KS, Omega, RhoCh, Bremsstrahlung, 10⁶ events per channel (486 + 327 s; Mac 190 + 118 s) | `spectra/New physics particles spectra/ALP-fermion/` |
| `13_nb3_sensitivity.sh` | nb3 SHiP-ECN3 (N≥2.3; 1271 s, Mac 382 s), FASER2-FPF (N≥2.3 and 3; 908 s, 316 s), FASER2-FPF-Etot100 (886 s, 279 s); production "Revised (2501.04525)", all 9 brem × DY variants | `Tabulated Nevents/`, `Sensitivity domains/ALP-fermion/<exp>/` |
| `14_curves_validate.sh` | curve files + validation (Python only, a few min): edges vs campaign/tracked/Maksym, deterministic channels (bremsstrahlung 0.5 %, fragmentation 2 %) per experiment, acceptance metadata, channel shares | `work/senscalc/out/`, `work/senscalc/validation/` |
| `15_noise_faser2.sh` | run-to-run Monte Carlo spread of FASER2-FPF: copy of the run tree, nb2 LHC/B only (10⁶) + nb3 FASER2-FPF | `work/senscalc/noise3/`, `validation/noise_faser2.txt` |

### Outputs of step 14 and where they would go

| File in `work/senscalc/out/` | Would replace | Check |
|---|---|---|
| `projections/faser2_senscalc_fpf.dat` (= `faser2_fpf_senscalc_nev3.dat`, AP-Central + DY-Central, N≥3) | `curves/data/bc10/projections/faser2_senscalc_fpf.dat` | edges vs tracked, 1 % tolerance (MC: not byte-identical) |
| `faser2_fpf_senscalc_nev2p3.dat` | campaign `out/faser2_fpf_senscalc_nev2p3.dat` (not tracked) | idem |
| `projections/ship_senscalc.dat` | `curves/data/bc10/projections/ship_senscalc.dat` | regenerated from M. Ovchynnikov's JSON (sha256 `70d6b4fc…dc2d` = the `curves/sources` file, reference tree only) with `curves.tools.senscalc_bc10_ship` itself: **byte-identical** (`97df2c6c…dded`); without the reference tree, a copy of the tracked file |
| `projections/ship_senscalc_rerun.dat` | nothing (validation of the pipeline) | vs Maksym, 0.5 % |

Validation targets, from the stored campaign outputs (reference tree; `py/compare_curves.py`, edges = crossings of the closed domains at
17 masses 0.05–4.5 GeV, log-log interpolation as `contour_tools.edges`):
- the campaign's own rerun vs Maksym's JSON, recomputed here from the stored files: lower edge +0.48 % at 3.5 GeV, upper edge
  −0.30 % at 4.0 GeV, mass range 0.0244–4.760 GeV identical in all 9 variants (the campaign README's "within 0.5 %");
- our rerun: SHiP vs Maksym (0.5 %) and vs the campaign rerun; FASER2-FPF and FASER2-FPF-Etot100 at N≥3 and N≥2.3, all
  9 variants, vs the campaign JSON (1 %);
- deterministic parts of the acceptance tables (Br_vis(m), cτ(m, y=1), facility numbers, grids) vs the campaign
  `logs/inspect_acceptances.txt`: must be identical (`py/compare_inspect.py`);
- channel shares at the FASER2 lower edge vs the campaign `logs/channel_breakdown_FASER2-FPF.txt`.

## Decoded tables (steps 20–23)

| Script | Inputs → outputs | Validation |
|---|---|---|
| `20_export_senscalc_2501.sh` | `third_party/SensCalc-v1.4` MX → `work/senscalc/senscalc_2501/` via the repository tool (2 min), plus `.gz` copies written by `py/gzip_zlib.py` (zlib level 9, mtime 0, OS 3 = macOS `gzip -n -9`; GNU gzip gives other bytes) | **done: all 5 tracked files byte-identical**, `.gz` included (`branching_ratios.csv.gz`, `decay_channels.json.gz`, `matrix_elements.json.gz`, `widths_bnt.csv.gz`, `widths_metadata.json`); would replace `grendel/models/alp/data/senscalc_2501` unchanged |
| `21_export_senscalc_2501_production.sh` | v1.4 MX → `work/senscalc/senscalc_2501_production/` | done; tool's hash/manifest checks pass. Against the only stored export (macOS 2026-07-13, SensCalc 0bca050, older exporter): production-coupling and fragmentation coefficients within 3e-14 (same inputs); meson files differ (the light-meson input changed between 0bca050 and v1.4); the stored fragmentation grid file is empty |
| `22_check_quartic.sh` | tracked `senscalc_quartic/*.dat` vs `phenomenology/Scalar/Production probabilities/` in eb2c456 and v1.4 | **done: all three tables and the LICENSE byte-identical** (no Mathematica step exists: SensCalc only reads these files) |
| `23_check_scalar_1809.sh` | decodes `phenomenology/Scalar/decay widths/ctauS.mx` (4 descriptions; the first, "1809.01876", is SensCalc's default) and the structure of `BrRatioScalar.mx` | done: tracked cτ (exHad copy) vs SensCalc "1809.01876" on 5508 common nodes: median 1.0e-4, 99 % ≤ 0.41 %, ≤ 1e-15 below 0.1 GeV, ≤ 0.83 % at 0.56–0.58 and 0.99 GeV, one outlier +38 % at the ππ-threshold node 0.280 GeV. Same tables up to small revisions between EventCalc/exHad and SensCalc: a provenance note, not a reproduction failure |

## Licence and parallelism

The free Wolfram Engine for developers needs a one-time activation with a Wolfram ID (done 2026-11-06, licence
`L9927-0674`, in `~/.WolframEngine/Licensing`). Limits: `$MaxLicenseProcesses` = 2, `$MaxLicenseSubprocesses` = 8, as for the
reference runs on macOS. The pipeline uses one main kernel + 4 subkernels (4 slots), never two SensCalc runs at once.
`$ProcessorCount` is 1 in this environment (12 cores), so the `Parallelization -> True` compiled functions of SensCalc run
single-threaded in each kernel: slower than on macOS, but each kernel stays on one core, within its slot. Results do not
depend on it (element-wise evaluation).

Slots: every step goes through `pipelines/common/slot_run.py -n K` (first come, first served since 2026-11-06 21:53):
K = 4 for a SensCalc run (main kernel + 4 subkernels), 1 for the light steps. `SENSCALC_SPARE_CORE=1` lets a 1-core step
under ~5 min (engine check, MX decoding) run on the core the 11-slot pool leaves free; nothing long runs outside the
pool. (The first production run, before the queue was first come, first served, used a 4-ticket workaround; it is gone.)
Memory: the campaign's main kernel peaked at 4.3 GB; `lib.sh` waits for 5 GB of available memory before each SensCalc run.

## Run

```
cd $REPRO && source env/activate.sh
bash pipelines/senscalc/00_install_wolfram.sh && bash pipelines/senscalc/01_fetch_sources.sh
wolframscript -activate                    # once, interactive, the user's Wolfram ID
setsid nohup bash pipelines/senscalc/run_all.sh > logs/senscalc/run_all.log 2>&1 < /dev/null &
```
Logs: `logs/senscalc/`. Summary of the validation: `work/senscalc/validation/summary.txt`.

## Results (2026-11-07)

Whole chain on 4 slots: about 1 h 45 min of running (nb1 40 min, nb2 14 min, nb3 51 min), 2.5–3.3× the macOS times
(single-threaded compiled code here), plus queueing in the shared pool. Summary: `work/senscalc/validation/summary.txt`.

| Product | Fig. | Result |
|---|---|---|
| `senscalc_2501` | 7 | **byte-identical**, all 5 tracked files including `.gz` |
| `senscalc_quartic` | 5 | **byte-identical** (upstream data files, no Mathematica step) |
| `ship_senscalc.dat` | 7 | regenerated from M. Ovchynnikov's JSON with the repository converter: **byte-identical** |
| SHiP-ECN3 rerun (pipeline validation) | – | vs campaign rerun ≤ 0.27 % in both edges at 17 masses, same mass range: **pass**; vs M. Ovchynnikov +0.56 % / −0.23 % (campaign: +0.48 % / −0.30 %) |
| `faser2_senscalc_fpf.dat` (`work/senscalc/out/projections/`) | 7 | not byte-identical (unseeded Monte Carlo); **within the measured run-to-run spread**, see below. Production run vs campaign: lower edge −2.13 % … +0.34 %, upper −2.07 % … +1.55 %; reach 3.717 vs 3.713 GeV. The 1 % criterion of step 14 is reported as FAIL for FASER2: a single 10⁶-event realisation (ours or the campaign's) does not meet it |
| Bremsstrahlung, fragmentation N_ev (step 14, every (m, y) with N ≥ 1) | – | FASER2-FPF and Etot100: bremsstrahlung ≤ 0.02 %, fragmentation ≤ 0.54 %: pass. SHiP-ECN3: bremsstrahlung ≤ 0.39 % (pass); fragmentation median 0.26 %, p99 1.9 %, max 2.6 % at one N ≈ 1 tail point (m = 2.3 GeV), reported as FAIL against the 2 % bar: SHiP's selection (E > 1 GeV, IP < 2.5 m) makes its acceptance partly Monte Carlo (the deterministic bremsstrahlung spectrum shows 0.3 % there against 0.002 % at FASER2 without cuts). No effect on the SHiP curve (≤ 0.27 %) |
| Acceptance metadata (Br_vis, cτ(m, y=1), facility numbers, grids) | – | identical to the campaign |
| scalar cτ cross-check (Figs. 4, 5) | – | tracked table vs SensCalc `ctauS.mx`: median 1e-4, one +38 % node at 0.280 GeV (provenance note) |

**FASER2-FPF run-to-run spread** (`15_noise_faser2.sh`, `validation/noise_faser2.txt`). Only the B → K(*) a spectrum is
Monte Carlo at the FASER2 edges below ~1 GeV (100 % of the lower edge at 0.1–0.5 GeV, 70–92 % of the upper edge at
0.3–0.5 GeV), drawn unseeded with 10⁶ events per channel; FASER2 (on axis, θ ≲ 1.5–2.3 mrad) sees only 0.15–0.4 % of
that sample, and its long-lifetime edge mostly the E < 300 GeV part (0.02–0.04 %). SHiP integrates a broad cone
(θ < 0.05: ours/campaign 0.996–0.998), hence its 0.27 %. Three independent B samples on Linux (run 1 = the pre-fix run,
whose B channel and lower edge are valid; run 2 = production; run 3 = noise run), N ≥ 3, central variant, 15 masses:
- lower edge: spread (max − min) median 1.0 %, max 3.1 % (0.2 GeV), sd 0.2–1.6 %; the campaign curve is within 0.78 %
  of the mean of the three at every mass, inside their min–max at 12/15 masses and ≤ 0.4 % outside at 0.3, 1.0, 2.0 GeV;
- upper edge (two samples; run 1 had no bremsstrahlung): spread median 0.8 %, max 4.6 % (0.4 GeV); the campaign is
  within 1.23 % of the mean at every mass (largest: −1.23 % at 3.5 GeV, −1.07 % at 0.7 GeV), inside the two-run range at
  4/11 masses; a two-sample range underestimates the spread.
The published curve is therefore one realisation of the same ~1 % noise; the offsets are not a systematic difference.

Effect on Fig. 7: the BC10 panel spans 6 decades in coupling over ~1000 px (220 dpi); 2 % is 0.0086 decade ≈ 1 px,
below the line width, and the largest lower-edge offset (0.2 GeV) lies left of the frame (x starts at 0.22 GeV). The
FASER2 curve does not visibly move.

**Smoke run** (2·10⁴ events): with 50× fewer events FASER2's forward B tail holds a few sampled ALPs per mass, and the
B-channel N_ev at the edges fluctuates by factors 0.004–15 against the campaign, so a curve tolerance of 5 % (the first
criterion) failed for statistical reasons alone (+32.7 % lower edge at 0.3 GeV, −20.4 % upper at 1.2 GeV). The
criterion now tests what a smoke run can: the chain completes, the deterministic channels are exact (bremsstrahlung
within 0.5 %: observed ≤ 0.029 %; fragmentation within 2 %: 0.69 %; this is the check the case-sensitivity bug fails,
at 0.001–0.03×), and the curve lands in the right place (edges within 50 %: +32.7 %; reach within 5 %: −1.3 %). Result:
**PASS** (`logs/senscalc/10_smoke_compare.log`).

Kept for the record: the pre-fix outputs in `work/senscalc/diag_case_bug_20261107/` (bremsstrahlung production
probabilities, JSONs, curve files, validation summary).
