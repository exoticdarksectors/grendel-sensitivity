# grendel-sensitivity

The sensitivity calculation behind **Proposing GRENDEL: GalleRy ExperimeNt
for Decays of Exotic LLPs at the LHC** (M. Citron, F. L. Redi, R. Schmitz,
[arXiv:2609.00152](https://arxiv.org/abs/2609.00152)): the machinery that
takes a Physics Beyond Colliders benchmark model from heavy-flavour
production at the HL-LHC to the exclusion island in its (mass, coupling)
plane, and draws that island over the published landscape.

Four benchmarks are implemented: the heavy neutral lepton with a single
dominant mixing (BC6–BC8, `|U|²`), the Higgs-portal dark scalar with mixing
only (BC4, `sin²θ`), the same scalar with the quartic coupling fixed by
BR(h→SS) = 0.01 (BC5, `sin²θ` from 0.14 GeV to m_h/2) and the fermiophilic
axion-like particle (BC10, `1/f`). They share one detector geometry, one
reconstruction, one acceptance Monte Carlo, one scan driver and one
uncertainty-band library; each model contributes its production, decay
model and numerical policy.

## Built on other people's work

This code is a thin layer over published tools and data, all fetched from
their own sources and named here so that credit lands where it belongs.
Cite them alongside the paper wherever they enter a result
(`THIRD_PARTY.md` has the pinned revisions, licences and every local
modification):

- **HNLCalc** — J. L. Feng, A. Hewitt, F. Kling, D. La Rocco,
  [arXiv:2405.07330](https://arxiv.org/abs/2405.07330): every HNL
  production branching ratio and three-body rate in this repository
  (mesons, baryons, tau) is computed by HNLCalc. We run a fork carrying
  six small fixes (proposed upstream), listed in `THIRD_PARTY.md`.
- **exHad** — V. Kryshtal, M. Ovchynnikov,
  [arXiv:2609.16104](https://arxiv.org/abs/2609.16104): the hadronic
  decay model behind the HNL, BC4 and BC10 rest-frame templates (with
  Pythia 8.3 and, for the HNL rates, [arXiv:1805.08567](https://arxiv.org/abs/1805.08567)).
- **FairShip** (SHiP Collaboration,
  [10.5281/zenodo.18020628](https://doi.org/10.5281/zenodo.18020628)) — the
  HNL width and decay-channel modules that produced the decay templates of
  the first published curve and the width-uncertainty table.
- **FONLL** — M. Cacciari, M. Greco, P. Nason
  ([hep-ph/9803400](https://arxiv.org/abs/hep-ph/9803400)); M. Cacciari,
  S. Frixione, P. Nason ([hep-ph/0102134](https://arxiv.org/abs/hep-ph/0102134)) with the
  **NNPDF4.0** parton distributions ([arXiv:2109.02653](https://arxiv.org/abs/2109.02653))
  and LHAPDF: the heavy-flavour production grids and their variations.
- **MadGraph5_aMC@NLO** — J. Alwall et al.,
  [arXiv:1405.0301](https://arxiv.org/abs/1405.0301): W/Z → ℓN with the **HeavyN**
  model of R. Ruiz and collaborators ([arXiv:1411.7305](https://arxiv.org/abs/1411.7305),
  [arXiv:1602.06957](https://arxiv.org/abs/1602.06957),
  [arXiv:1812.08750](https://arxiv.org/abs/1812.08750)), and the SM
  prompt-tau pool.
- **Pythia 8** — C. Bierlich et al.,
  [arXiv:2203.11601](https://arxiv.org/abs/2203.11601): the kaon spectrum
  of the K → ℓN channel, the BC10 decay templates (`templates_pythia`) and
  the event generation inside exHad and FairShip.
- **SensCalc** — M. Ovchynnikov et al.,
  [arXiv:2305.13383](https://arxiv.org/abs/2305.13383), with the
  fermiophilic-ALP phenomenology of [arXiv:2501.04525](https://arxiv.org/abs/2501.04525):
  the BC10 widths and branching ratios.
- **ALPINIST** — J. Jerhot, B. Döbrich, F. Ertas, F. Kahlhoefer, T. Spadaro,
  [arXiv:2201.05170](https://arxiv.org/abs/2201.05170): the BC10 production
  (its RG evaluation of the b → s a coupling of GKOZ,
  [arXiv:2310.03524](https://arxiv.org/abs/2310.03524), ported in
  `grendel/models/alp/tools/compute_cbs_alpinist.py`, and its B → K⁽ⁱ⁾ a
  form factors), the NA62 bound and the digitised GKOZ widths (audit only).
- **HNLimits** — E. Fernández-Martínez, M. González-López, J. Hernández-García,
  M. Hostert, J. López-Pavón, [arXiv:2304.06772](https://arxiv.org/abs/2304.06772):
  the experimental bounds on the HNL panels.
- **Y.-M. Chen, Y. Zhang**, [arXiv:2410.07343](https://arxiv.org/abs/2410.07343)
  (Phys. Rev. D 111 (2025) 123024): the BBN bound on the HNL panels, computed
  for a single Majorana HNL (95 % CL).
- **A. Escalante del Valle**, [arXiv:2509.02564](https://arxiv.org/abs/2509.02564)
  (JHEP 12 (2025) 129): the ATLAS/CMS h → S S long-lived-particle limits
  in the BC5 plane, drawn from its published data.
- **M. W. Winkler**, [arXiv:1809.01876](https://arxiv.org/abs/1809.01876):
  the dark-scalar hadronic widths and the inclusive b → X_s S rate; above
  7.5 GeV the same model continued by the EventCalc tables shipped with exHad
  (`grendel/models/scalar/data/eventcalc_1809`).
- **Boiarska et al.**, [arXiv:1904.10447](https://arxiv.org/abs/1904.10447)
  and SensCalc: the quartic-coupling production of BC5 (h → SS, B → X_s SS,
  B_s → SS; tables in `grendel/models/scalar/data/senscalc_quartic`).
- The **LHC Higgs Cross Section Working Group** (YR4,
  [arXiv:1610.07922](https://arxiv.org/abs/1610.07922)): the 14 TeV Higgs
  cross sections that normalise the BC5 h → SS channel, whose kinematics
  come from MadGraph5_aMC@NLO (LO, every SM mode) showered with Pythia 8.
- The **PBC** reports ([arXiv:1901.09966](https://arxiv.org/abs/1901.09966),
  [arXiv:2310.17726](https://arxiv.org/abs/2310.17726),
  [arXiv:2505.00947](https://arxiv.org/abs/2505.00947)) and the SHiP, FASER2,
  CODEX-b, ANUBIS and DarkQuest papers whose projections are digitised in
  `curves/data`, each named in its file's header or its directory's `SOURCE.md`.

## Layout

```
grendel/
  constants.py   luminosity, the N_sig >= 3 criterion, the HNL |U|² grid
  geometry/      the detector mesh, the ray-cast cache and the two-track vertex reconstruction
  reco/          the acceptance Monte Carlo, the exclusion-island solver, decay-template tooling
  scan.py        the per-mass scan driver every model runs through (ModelSpec contract, resume, thresholds)
  cutflow.py     the selection cutflow at a chosen coupling, from the same Monte Carlo
  production/    FONLL heavy-flavour grids, the meson sampler and the FONLL grid generator
  models/        hnl, scalar (BC4 and BC5), alp (BC10): production, decay templates, ModelSpec, plots
  band/          the uncertainty campaigns (FONLL scale/PDF/mass, decay model, controls) and their combination
  io/            paths policy (environment variables), file formats, durable writes
curves/          reference curves of other experiments, the figure renderer, the digitisation tools
third_party/     fetcher for the external software, pinned (exHad, HNLCalc, FairShip, HeavyN UFO, FONLL, MadGraph, Pythia)
scripts/         consistency checks
tests/           unit tests and an end-to-end smoke pipeline
```

## Tunnel geometry (since October 2026)

The tunnel is UP542 **as built** (measured section, as-built centreline and floor, official CMS frame, correct
handedness), mirrored unchanged from the single geometry source https://gitlab.cern.ch/grendel/grendel-geometry in
`grendel/geometry/data/up542/`; see [`GEOMETRY_UP542.md`](GEOMETRY_UP542.md). The tunnel of the published curves
(v1–v2.2) is in those tags; from v2.3 on every command below uses the UP542 tunnel of grendel-geometry v1.0. The geometry files are
CERN-internal: this repository is private.

## Install

```
git clone https://github.com/exoticdarksectors/grendel-sensitivity
cd grendel-sensitivity
python -m venv .venv && . .venv/bin/activate
pip install -e ".[curves,test]"        # Python >= 3.11
python -m pytest tests/unit -q          # ~10 s
python -m pytest tests/e2e -q           # the whole chain at tiny statistics, ~15 s
```

`matplotlib` is pinned to 3.10.9 because the paper figures reproduce byte
for byte with it; `HNLimits` (the `curves` extra) provides the HNL bounds.
The optional extras are `raycast` (embreex, the accelerated ray backend the
campaigns used), `production` (the HNL production dependencies) and
`digitize` (PyMuPDF, for re-extracting reference curves from PDFs).

## From production to figure

Every stage is a `python -m` entry point with `--help`. Inputs and outputs
resolve through `grendel.io.paths`: set `GRENDEL_WORK_DIR` (default
`./grendel_work`) or the per-model `GRENDEL_<HNL|BC4|BC5|BC10>_{VECTORS,TEMPLATES,ANALYSIS,GEOMETRY}_DIR`,
or pass the directories on the command line.

**1. External software.** `python third_party/fetch.py --list` names what
each stage needs and where it comes from; `--all` fetches everything at the
pinned revisions. BC4 and BC10 run on the tracked FONLL grids alone; the HNL
production needs HNLCalc (and MadGraph + LHAPDF for the W/Z and prompt-tau
channels); every published decay-template generator needs either exHad in
its own virtual environment or PyROOT + Pythia 8 (the legacy BC10 proxy
needs neither).

**2. Production** — four-vectors per mass point.

```
python -m grendel.models.scalar.scan --help          # BC4 produces its own (FONLL bottom -> b-hadron -> X_s S, inclusive)
python -m grendel.models.scalar.scan_bc5 --help      # BC5: the BC4 rows plus h -> S S and b -> X_s S S
python -m grendel.models.alp.production               # BC10: B -> K a, K* a ... from the FONLL grid
python -m grendel.models.hnl.production.run_all       # HNL: D, B, Bc, Lambda_b, tau, K and W/Z channels, combined
```

**3. Decay templates** — rest-frame final states per mass point.

```
python -m grendel.models.hnl.templates.exhad --flavor Ue Umu Utau      # exHad (needs GRENDEL_EXHAD_PYTHON's venv)
python -m grendel.models.hnl.templates.fairship                        # the FairShip model of the first curve (PyROOT)
python -m grendel.models.scalar.templates_exhad --model scalar-1809
python -m grendel.models.alp.templates_exhad | templates_pythia | templates_proxy --legacy-proxy
```

BC4 also has an analytic two-body decay backend that needs no templates, and
BC5 falls back to it over the whole range without `--templates-dir` (the
b bbar continuum above 2 m_B then enters as a two-track proxy like s s, c c
and g g). The generators of one model are alternatives writing the same
`GRENDEL_<model>_TEMPLATES_DIR` (pass `--out` to keep several); the HNL and
BC10 scans read that directory, the BC4 and BC5 scans only with
`--templates-dir`. The v2 curves use exHad for all four benchmarks; the
exact commands (the scalar set over the BC5 grid, BC10 in two passes) are in
`curves/data/grendel_v2/README.md`.

The BC5 Higgs kinematics are a tracked histogram
(`grendel/models/scalar/data/higgs_pt_y_14TeV.csv.gz`) built once by three
stages that need MadGraph + LHAPDF and a Python with `pythia8`:
`grendel.models.scalar.higgs_production` (LO hard process of ggF, VBF, WH,
ZH, ttH), `higgs_shower` (Pythia 8 parton shower) and `higgs_pool` (the
histogram, normalised per mode to the YR4 14 TeV cross sections). The scan
itself reads only the histogram.

**4. Scan** — the exclusion island, solved at one or several thresholds.

```
python -m grendel.models.hnl.scan    --flavor Ue Umu Utau --templates-dir <exHad HNL templates> --thresholds 3 10 --out RESULTS/hnl
python -m grendel.models.scalar.scan --n-pool 400000 --high-pt-tilt-scale 5 --templates-dir <exHad scalar templates> --thresholds 3 10 --out RESULTS/bc4
python -m grendel.models.scalar.scan_bc5 --templates-dir <exHad scalar templates> --thresholds 3 10 --out RESULTS/bc5
python -m grendel.models.alp.scan    --templates-dir <exHad BC10 templates> --thresholds 3 10 --out RESULTS/bc10
python -m grendel.io.thresholds RESULTS/<model>/sensitivity.csv --threshold 10   # -> sensitivity_nsig10.csv
python -m grendel.io.variants RESULTS/bc5/sensitivity.csv --variant hSS          # BC5 by mode; --list names them all
```

These are the v2 settings (`--n-pool 400000 --high-pt-tilt-scale 5` is the
published BC4 pool, not the CLI default); without `--templates-dir` BC4 and
BC5 use the analytic decays of v1.

Each run writes `sensitivity.csv`, `run_metadata.json`, `scan_status.json`
and a `geometry_cache/` it reuses; `--resume` keeps the rows computed with
the same settings, code and input files (SHA-256 recorded per row). The
selection cutflow at any point of the island is
`python -m grendel.cutflow --model bc4 --mass 1.0 --from-results RESULTS/bc4/sensitivity.csv`
(`--model bc5` likewise).

*BC5 by production mode (v2.1).* BC5 makes the scalar three ways: `mixing`
(b → X_s S through the mixing angle, the BC4 channel), `hSS` (on-shell
h → SS, which only an LHC experiment has) and `BSS` (b → X_s SS and
B_s → SS through the off-shell Higgs, the only Higgs-mediated piece a beam
dump sees). Besides the nominal island, every row of the BC5
`sensitivity.csv` carries the island of each mode alone (`mixing_u2_min`,
`hSS_peak_N`, `BSS_has_sensitivity`, ...) and, for each `--br-hss-overlay`
value (default 0.001), the islands rescaled to that BR(h → SS): the full
curve (`brhss0.001_*`) and each quartic mode (`brhss0.001_hSS_*`,
`brhss0.001_BSS_*`). The quartic yields are linear in BR(h → SS), so all of
these come from the same Monte Carlo as the nominal curve.
`grendel.io.variants` writes any of them in the layout of `sensitivity.csv`,
and `grendel.io.thresholds` applies to it as to any table;
`scripts/plot_bc5_modes.py` draws the split and the BR scaling. The mode of
every four-vector row is the seventh column of the BC5 CSVs; CSVs written by
v2 lack it, and the scan regenerates them (same seeds, same pools).

**5. Uncertainty bands** (optional; needs the FONLL variation campaign of
`grendel.production.fonll_grids`): `grendel.band.scalar.campaign`,
`grendel.band.alp.campaign` and `grendel.band.hnl.campaign` re-run the scan
over the production variations; `campaign collect` (BC4, BC10) or
`grendel.band.quadrature` (HNL) builds the band;
`grendel.band.hnl.{decay_model_band,bc_nuisance,decay_model_band_exhad}`
add the HNL decay-model and B_c nuisances. See `grendel/band/README.md`.

**6. Figures.**

```
python -m curves.plot paper --grendel-dir RESULTS --out figures
```

reads `RESULTS/<model>/sensitivity{,_nsig10}.csv` and writes the manuscript
figures (HNL, BC4, BC5, BC10); `talk` and `diagnostic` are the other layouts.
The frozen tables of each release render the same way:
`--grendel-dir curves/data/grendel_v1` (arXiv:2609.00152v1) and
`curves/data/grendel_v2` (v2). The figures of the second version of the paper
(arXiv:2609.00152v2), on the UP542 tunnel, come from

```
python -m curves.plot paper --grendel-dir curves/data/grendel_v2_5 --benchmark bc4 bc5
python -m curves.plot paper --grendel-dir curves/data/grendel_v2_4 --benchmark bc10
python -m curves.plot paper --grendel-dir curves/data/grendel_v2_3 --benchmark hnl     # HNLimits==1.2.0
```

(tag v2.5.0; BC4 and BC5 with the exHad/EventCalc scalar model of the SHiP
projections, the HNL with a fine grid at 2.0-4.2 GeV, BC10 with
the arXiv:2310.03524 coupling and the ALP widths of exHad e41b07b; the READMEs
of `grendel_v2_3`, `grendel_v2_4` and `grendel_v2_5` have the recipes, hashes and checks). The same panels on the
arXiv:2609.00152 tunnel come from `grendel_v2_1_hybrid` (BC4, BC5),
`grendel_v2_refined` (HNL) and `grendel_v2_2` (BC10), tag v2.2.2. The
comparison curves, their sources and the conventions of every panel are
documented in `curves/README.md` and `curves/PROVENANCE.md`.

## What the numbers depend on

The published curves are fixed by choices this repository states rather
than hides: the ray-cast geometry cache, the per-mass random seeds
(`grendel/models/{hnl,scalar,alp}/spec.py`), the event-chunking policy, the exact-
versus-resampled hit selection, the exclusion solver (`grendel/reco/exclusion.py`)
and the decay model of the templates. `run_metadata.json` records the
settings of every scan. The published curves themselves are not part of
this repository. The published HNL `N_sig >= 3` curve is a splice of two
historical runs and is not bit-reproducible from one checkout: a
self-consistent re-run differs on 40 of 369 rows (median 0.11 %, max 2.7 %).

Not everything can be re-run from a clean checkout: the HNL production
needs HNLCalc, MadGraph and LHAPDF; the template generators need exHad or
PyROOT; the FONLL variation campaign needs a FONLL build with LHAPDF and
days of CPU; the campaigns' accelerated ray backend needs embreex. Each of
these is fetched or installed by the reader, never redistributed here, and
`THIRD_PARTY.md` lists the terms and the modifications the pinned revisions
carry. The FONLL modifications ship as diffs
(`grendel/production/fonll_grids/patches`).

## Citing

Cite the paper (`CITATION.cff`, which also lists the references of the
software above) and, for the inputs that contribute to a result, the
packages named under "Built on other people's work" and the sources in the
headers of the reference curves.

## Licence

BSD-3-Clause (`LICENSE`; `grendel/models/alp/tools/compute_cbs_alpinist.py`
is a port of ALPINIST code under ALPINIST's BSD-3-Clause,
`curves/sources/alpinist/LICENSE`), except two files adapted from FairShip —
`grendel/models/hnl/templates/fairship_decay.py` and
`grendel/models/hnl/data/fairship_decay_selection.conf` — which carry
FairShip's LGPL-3.0-or-later (`LICENSES/LGPL-3.0-or-later.txt`); each states
so in its own header. The external software the pipeline drives keeps its
own licences; see `THIRD_PARTY.md`.
