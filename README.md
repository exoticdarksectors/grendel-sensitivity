# grendel-sensitivity

The sensitivity calculation behind **Proposing GRENDEL: GalleRy ExperimeNt
for Decays of Exotic LLPs at the LHC** (M. Citron, F. L. Redi, R. Schmitz,
[arXiv:2609.00152](https://arxiv.org/abs/2609.00152)): the machinery that
takes a Physics Beyond Colliders benchmark model from heavy-flavour
production at the HL-LHC to the exclusion island in its (mass, coupling)
plane, and draws that island over the published landscape.

Three benchmarks are implemented: the heavy neutral lepton with a single
dominant mixing (BC6–BC8, `|U|²`), the Higgs-portal dark scalar (BC4,
`sin²θ`) and the fermiophilic axion-like particle (BC10, `1/f`). They share
one detector geometry, one reconstruction, one acceptance Monte Carlo,
one scan driver and one uncertainty-band library; each model contributes
its production, decay model and numerical policy.

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
  four small fixes, listed in `THIRD_PARTY.md`.
- **exHad** — V. Kryshtal, M. Ovchynnikov,
  [arXiv:2609.16104](https://arxiv.org/abs/2609.16104): the hadronic
  decay model behind the HNL, BC4 and BC10 rest-frame templates (with
  Pythia 8.3 and, for the HNL rates, [arXiv:1805.08567](https://arxiv.org/abs/1805.08567)).
- **FairShip** (SHiP Collaboration) — the HNL width and decay-channel
  modules that produced the decay templates of the first published curve
  and the width-uncertainty table.
- **FONLL** — M. Cacciari, S. Frixione, P. Nason
  ([hep-ph/9803400](https://arxiv.org/abs/hep-ph/9803400),
  [hep-ph/0102134](https://arxiv.org/abs/hep-ph/0102134)) with the
  **NNPDF4.0** parton distributions ([arXiv:2109.02653](https://arxiv.org/abs/2109.02653))
  and LHAPDF: the heavy-flavour production grids and their variations.
- **MadGraph5_aMC@NLO** — J. Alwall et al.,
  [arXiv:1405.0301](https://arxiv.org/abs/1405.0301), with the **HeavyN**
  model of R. Ruiz and collaborators ([arXiv:1411.7305](https://arxiv.org/abs/1411.7305),
  [arXiv:1602.06957](https://arxiv.org/abs/1602.06957),
  [arXiv:1812.08750](https://arxiv.org/abs/1812.08750)): W/Z and
  prompt-tau production of HNLs.
- **Pythia 8** — C. Bierlich et al.,
  [arXiv:2203.11601](https://arxiv.org/abs/2203.11601): the kaon spectrum
  of the K → ℓN channel and the event generation inside exHad and FairShip.
- **SensCalc** — M. Ovchynnikov et al.,
  [arXiv:2305.13383](https://arxiv.org/abs/2305.13383), with the
  fermiophilic-ALP phenomenology of [arXiv:2501.04525](https://arxiv.org/abs/2501.04525)
  and [arXiv:2310.03524](https://arxiv.org/abs/2310.03524) (GKOZ): the BC10
  widths and branching ratios.
- **ALPINIST** — J. Jerhot, B. Döbrich, F. Ertas, F. Kahlhoefer, T. Spadaro,
  [arXiv:2201.05170](https://arxiv.org/abs/2201.05170): the NA62 bound and
  the digitised GKOZ widths.
- **HNLimits** — E. Fernández-Martínez, M. González-López, J. Hernández-García,
  M. Hostert, J. López-Pavón, [arXiv:2304.06772](https://arxiv.org/abs/2304.06772):
  the experimental and cosmological bounds on the HNL panels.
- **M. W. Winkler**, [arXiv:1809.01876](https://arxiv.org/abs/1809.01876):
  the dark-scalar hadronic widths.
- The **PBC** reports ([arXiv:1901.09966](https://arxiv.org/abs/1901.09966),
  [arXiv:2310.17726](https://arxiv.org/abs/2310.17726),
  [arXiv:2505.00947](https://arxiv.org/abs/2505.00947)) and the SHiP, FASER2,
  CODEX-b, ANUBIS and DarkQuest papers whose projections are digitised in
  `curves/data`, each named in its file's header.

## Layout

```
grendel/
  constants.py   luminosity, the N_sig >= 3 criterion, the coupling grid
  geometry/      the detector mesh, the ray-cast cache and the two-track vertex reconstruction
  reco/          the acceptance Monte Carlo, the exclusion-island solver, decay-template tooling
  scan.py        the per-mass scan driver every model runs through (ModelSpec contract, resume, thresholds)
  cutflow.py     the selection cutflow at a chosen coupling, from the same Monte Carlo
  production/    FONLL heavy-flavour grids, the meson sampler and the FONLL grid generator
  models/        hnl, scalar (BC4), alp (BC10): production, decay templates, ModelSpec, plots
  band/          the uncertainty campaigns (FONLL scale/PDF/mass, decay model, controls) and their combination
  io/            paths policy (environment variables), file formats, durable writes
curves/          reference curves of other experiments, the figure renderer, the digitisation tools
third_party/     fetcher for the external software, pinned (exHad, HNLCalc, FairShip, HeavyN UFO, FONLL, MadGraph, Pythia)
scripts/         consistency checks
tests/           unit tests and an end-to-end smoke pipeline
```

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
`./grendel_work`) or the per-model `GRENDEL_<HNL|BC4|BC10>_{VECTORS,TEMPLATES,ANALYSIS,GEOMETRY}_DIR`,
or pass the directories on the command line.

**1. External software.** `python third_party/fetch.py --list` names what
each stage needs and where it comes from; `--all` fetches everything at the
pinned revisions. BC4 and BC10 run on the tracked FONLL grids alone; the HNL
production needs HNLCalc (and MadGraph + LHAPDF for the W/Z and tau
channels); every decay-template generator needs either exHad in its own
virtual environment or PyROOT + Pythia 8.

**2. Production** — four-vectors per mass point.

```
python -m grendel.models.scalar.scan --help          # BC4 produces its own (FONLL bottom -> B -> K S)
python -m grendel.models.alp.production               # BC10: B -> K a, K* a ... from the FONLL grid
python -m grendel.models.hnl.production.run_all       # HNL: D, B, Bc, Lambda_b, tau, K and W/Z channels, combined
```

**3. Decay templates** — rest-frame final states per mass point.

```
python -m grendel.models.hnl.templates.exhad --flavor Ue Umu Utau      # exHad (needs GRENDEL_EXHAD_PYTHON's venv)
python -m grendel.models.hnl.templates.fairship                        # the FairShip model of the first curve (PyROOT)
python -m grendel.models.scalar.templates_exhad --model scalar-1809
python -m grendel.models.alp.templates_exhad | templates_pythia | templates_proxy
```

BC4 also has an analytic two-body decay backend that needs no templates.

**4. Scan** — the exclusion island, solved at one or several thresholds.

```
python -m grendel.models.hnl.scan    --flavor Ue Umu Utau --thresholds 3 10 --out RESULTS/hnl
python -m grendel.models.scalar.scan --thresholds 3 10 --out RESULTS/bc4
python -m grendel.models.alp.scan    --thresholds 3 10 --out RESULTS/bc10
python -m grendel.io.thresholds RESULTS/<model>/sensitivity.csv --threshold 10   # -> sensitivity_nsig10.csv
```

Each run writes `sensitivity.csv`, `run_metadata.json`, `scan_status.json`
and a `geometry_cache/` it reuses; runs resume with `--resume`. The
selection cutflow at any point of the island is
`python -m grendel.cutflow --model bc4 --mass 1.0 --from-results RESULTS/bc4/sensitivity.csv`.

**5. Uncertainty bands** (optional; needs the FONLL variation campaign of
`grendel.production.fonll_grids`): `grendel.band.scalar.campaign`,
`grendel.band.alp.campaign` and `grendel.band.hnl.campaign` re-run the scan
over the production variations and `... combine` builds the envelope;
`grendel.band.hnl.{decay_model_band,bc_nuisance,decay_model_band_exhad}`
add the HNL decay-model and B_c nuisances. See `grendel/band/README.md`.

**6. Figures.**

```
python -m curves.plot paper --grendel-dir RESULTS --out figures
```

reads `RESULTS/<model>/sensitivity{,_nsig10}.csv` and writes the three
manuscript figures; `talk` and `diagnostic` are the other layouts. The
comparison curves, their sources and the conventions of every panel are
documented in `curves/README.md` and `curves/PROVENANCE.md`.

## What the numbers depend on

The published curves are fixed by choices this repository states rather
than hides: the ray-cast geometry cache, the per-mass random seeds
(`grendel/models/<model>/spec.py`), the event-chunking policy, the exact-
versus-resampled hit selection, the exclusion solver (`grendel/reco/exclusion.py`)
and the decay model of the templates. `run_metadata.json` records the
settings of every scan. The published curves themselves are not part of
this repository.

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

BSD-3-Clause (`LICENSE`). The external software the pipeline drives keeps
its own licences; see `THIRD_PARTY.md`.
