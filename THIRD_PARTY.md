# Third-party software and data

This repository redistributes no external software. The programs and
models the pipeline drives are fetched from their upstreams at pinned
revisions by `third_party/fetch.py` (see `third_party/README.md`), and each
keeps its own licence. This file records what is used, under which terms,
and every modification the pinned revisions carry relative to their
origins. Data derived from published work is listed at the end with the
citation each carries.

## Fetched software

### exHad — V. Kryshtal, M. Ovchynnikov — BSD-3-Clause

`maksymovchynnikov/exHad` at commit `7bb2308d1a60090603cbf07298d15b039f750727`
(2026-09-16). V. Kryshtal, M. Ovchynnikov, *Rethinking search signatures:
hadronic decays of GeV-scale feebly coupled particles*, arXiv:2609.16104.
Cite it, the rate calculations its README lists for each model (for the
HNL: arXiv:1805.08567) and Pythia 8.3 (arXiv:2203.11601) wherever exHad
decays are used. exHad's own `LICENSE`, `LICENSE-EventCalc` (the EventCalc
kinematics and rate tables it incorporates, BSD-3-Clause) and
`THIRD_PARTY.md` (DeLiVeR, ReD-DeLiVeR, the numerical inputs) apply to its
checkout. Used unmodified, in its own virtual environment, to generate the
rest-frame decay templates of the HNL, dark-scalar and fermiophilic-ALP
models (`grendel.reco.exhad_generation`).

### HNLCalc — J. L. Feng, A. Hewitt, F. Kling, D. La Rocco — no licence published

HNLCalc computes every HNL production branching ratio and three-body
rate used here; cite J. L. Feng, A. Hewitt, F. Kling, D. La Rocco,
*Simulating Heavy Neutral Leptons with General Couplings at Collider and
Fixed Target Experiments*, arXiv:2405.07330, wherever an HNL production
rate enters a result. Upstream `laroccod/HNLCalc` carries no licence file, which is why it is fetched and
never vendored. The pinned revision is the fork `leoredi/HNLCalc` at
`e292cef9c82964b493fda5d3b3be90a8fa44b235`, which is upstream `main` at
`07f84728` plus four commits:

1. the unused `from skhep.math.vectors import LorentzVector, Vector3D` is
   commented out (it forced a dependency on a module dropped from modern
   scikit-hep; no numerical effect);
2. a duplicated `Ds+ -> K0` form-factor block is removed: two plain `if`
   statements matched the same condition and assigned `f00 = 0.747` then
   `f00 = 0.72`, so the first was dead code; the effective `0.72` (Melikhov
   and Stech, hep-ph/0001113, Table XVII) is kept and the behaviour is
   unchanged;
3. a misplaced parenthesis in `get_2body_br_tau` (pseudoscalar branch,
   `tau -> P N`) is fixed: the Källén factor was `sqrt(1 - A*(1 - B))`
   instead of `sqrt((1 - A)*(1 - B))`, overestimating `BR(tau -> pi N)` by
   x1.07 at m_N = 0.5 GeV, x1.38 at 1.0 GeV and x3.9 at 1.5 GeV. The
   upstream form is still present in the public repository; this is a
   deliberate local deviation;
4. a double-squared mixing factor in `get_3body_dbr_baryon` is fixed
   (`|U|^4` -> `|U|^2`); numerically inert for the single-flavour instances
   used here, wrong for mixed couplings.

Used for the meson, b-baryon and tau production branching ratios and
three-body differential rates (`grendel.models.hnl.production`); the
lifetimes and visible final states come from the decay templates, not
from HNLCalc.

### FairShip — the SHiP Collaboration — LGPL-3.0-or-later

`ShipSoft/FairShip`, only its `python/` directory, at the revision in
`PINS.json` (March 2026). Copyright CERN for the benefit of the SHiP
Collaboration. The modules `hnl.py`, `readDecayTable.py`,
`pythia8_conf_utils.py`, `shipunit.py` and `alpha_s.dat` supply the
flavour-dependent HNL widths and the Pythia 8 decay-channel configuration
behind the FairShip decay templates of the first published curve and the
width-band table (`grendel.models.hnl.templates.fairship_decay`,
`grendel.band.hnl.width_band`). Used unmodified: the pipeline puts the
directory on `sys.path` and imports ROOT itself before loading them.
Needs PyROOT with Pythia 8 (`ROOT.TPythia8`); FairShip itself is not built.

### HeavyN UFO model `SM_HeavyN_CKM_AllMasses_LO` — R. Ruiz and collaborators — FeynRules model database

Version 2.2 (R. Ruiz, 27 June 2016), downloaded from the FeynRules HeavyN
model page, archive SHA-256 in `PINS.json`. C. Degrande, O. Mattelaer,
R. Ruiz, J. Turner, arXiv:1602.06957; S. Pascoli, R. Ruiz, C. Weiland,
arXiv:1812.08750. Used with MadGraph for W/Z -> l N production
(`grendel.models.hnl.production.madgraph.run_wz_production`). MadGraph 3
rewrites `object_library.py` and `write_param_card.py` for Python 3 when
it first imports the model; the model's physics files are untouched.

### FONLL 1.3.3 — M. Cacciari, S. Frixione, P. Nason — distributed by its authors

M. Cacciari, S. Frixione, P. Nason. Fetched from the `alisw/fonll` mirror at
`c7086e49141cf6705cf7a4bc5f7d0b3a38673203` (the 1.3.3 import, identical to
the authors' tarball). The modifications the grids were produced with ship
as unified diffs in `grendel/production/fonll_grids/patches` (rapidity
capacity, BCFY fragmentation modes, output precision, the LHAPDF link
target) and are applied by `grendel.production.fonll_grids.install`; see
that directory's README for each hunk. Cite Cacciari, Greco, Nason, JHEP
9805 (1998) 007 and Cacciari, Frixione, Nason, JHEP 0103 (2001) 006, as
the FONLL web form asks. Keep FONLL's attribution intact; obtain it
through the official channel.

### MadGraph5_aMC@NLO 3.6.6 — J. Alwall et al. — UoI-NCSA (MadGraph licence)

`mg5amcnlo/mg5amcnlo` tag `v3.6.6`. Alwall et al., JHEP 07 (2014) 079,
arXiv:1405.0301. Used for the W/Z -> l N and prompt-tau production
samples. Needs LHAPDF 6 with `NNPDF40_nlo_as_01180`.

### Pythia 8.315 and 8.317 — C. Bierlich et al. — GPL-2.0-or-later

`Pythia8/releases` tags `pythia8315` and `pythia8317`. Bierlich et al.,
SciPost Phys. Codebases 8 (2022), arXiv:2203.11601. 8.315 reproduces the
tracked kaon spectrum (`grendel.models.hnl.production.make_kaon_spectrum`;
8.317 shifts <n_K> at the 0.5 % level); 8.317 is the version exHad
requires.

### Python packages

Installed from PyPI by `pip` at the versions in `pyproject.toml`. The
figures need `HNLimits==1.2.0` (mhostert/Heavy-Neutrino-Limits, MIT;
E. Fernández-Martínez, M. González-López, J. Hernández-García, M. Hostert,
J. López-Pavón, arXiv:2304.06772), whose bundled data files are the
experimental bounds on the HNL panels — cite that paper with the panels;
`matplotlib==3.10.9` is pinned because the manuscript figures reproduce
byte for byte with it.

## Data derived from published work

- `grendel/models/alp/data/senscalc_2501/` — fermiophilic-ALP decay widths
  and branching ratios exported from SensCalc (M. Ovchynnikov, J.-L. Tastet,
  O. Mikulenko, K. Bondarenko, arXiv:2305.13383; BSD-3-Clause), the revised phenomenology of arXiv:2501.04525; its
  `PROVENANCE.md` records the export.
- `grendel/models/alp/data/alpinist/` — digitised GKOZ (arXiv:2310.03524)
  width curves from ALPINIST (J. Jerhot, B. Döbrich, F. Ertas, F. Kahlhoefer,
  T. Spadaro, arXiv:2201.05170; BSD-3-Clause; licence in
  `curves/sources/alpinist/LICENSE`), retained for audit only.
- `grendel/models/scalar/data/winkler_widths.csv` — dark-scalar hadronic
  widths digitised from M. W. Winkler, Phys. Rev. D 99 (2019) 015018,
  arXiv:1809.01876, Fig. 4.
- `curves/data/` — reference curves digitised from the figures of the
  papers named in each file's header (PBC reports arXiv:1901.09966,
  arXiv:2310.17726, arXiv:2505.00947; SHiP, FASER2, CODEX-b, ANUBIS,
  DarkQuest, GKOZ, FIPs 2022, Sabti et al.), the HNLimits workbook, and
  the ALPINIST NA62 bound; see `curves/README.md` and `curves/PROVENANCE.md`.
- `grendel/production/data/fonll/central/` — meson grids produced with the
  software above (`grendel/production/data/fonll/PROVENANCE.md`); cite
  FONLL and NNPDF4.0 (arXiv:2109.02653) when using them.
- `grendel/geometry/` — the detector geometry and reconstruction shared
  with the Higgs-portal study of the same paper; `grendel/geometry/PROVENANCE.md`.
