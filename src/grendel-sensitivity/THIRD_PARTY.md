# Third-party software and data

This repository ships no copy of an external program: the software the
pipeline drives is fetched from its upstream at a pinned revision by
`third_party/fetch.py` (see `third_party/README.md`), and each keeps its own
licence. The exceptions are two small FairShip-derived files, listed under
FairShip below, which carry FairShip's LGPL-3.0-or-later, and
`grendel/models/alp/tools/compute_cbs_alpinist.py`, a port of ALPINIST's RG
code under ALPINIST's BSD-3-Clause (`curves/sources/alpinist/LICENSE`). This file records
what is used, under which terms, and every modification the pinned revisions
carry relative to their origins. Data derived from published work is listed
at the end with the citation each carries.

## Fetched software

### exHad — V. Kryshtal, M. Ovchynnikov — BSD-3-Clause

`maksymovchynnikov/exHad` at commit `e41b07b2ba8004e7f271ac0e9920d47819131b99`
(2026-10-16; the ALP-fermion tables with the corrected renormalization-group running; the
scalar and HNL data are those of `7bb2308`, 2026-09-16, behind the earlier tables). V. Kryshtal, M. Ovchynnikov, *Rethinking search signatures:
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
`b1d36c6bc429f3041c730ef77ce8347c5124e694` (its `main`, proposed upstream as
laroccod/HNLCalc#6), which is upstream `main` at `07f84728` plus six fixes.
Up to v2.2 the pin was `e292cef` (fixes 1-4 only; kept by the fork's tag
`grendel-v2.2`).

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
   used here, wrong for mixed couplings;
5. one Fermi constant, the PDG value 1.1663788e-5 GeV^-2, replaces the
   per-method definitions (1.166378e-5, 1.1663787e-5) and the unused `GF()`
   (laroccod/HNLCalc#4); effect ~1e-6;
6. f_D and f_Ds are those of FLAG 2024 (arXiv:2411.04268, N_f = 2+1+1),
   212.0 and 249.9 MeV, instead of 222.6 and 280.1 MeV: D_s -> l N rates
   -20 %, D -> l N -9 %; BR(D_s -> tau nu) is now 5.24 % (BESIII 2023:
   5.44(21) %), and D_s -> l N agrees with SensCalc v1.4 to 2 %.

Used for the meson, b-baryon and tau production branching ratios and
three-body differential rates (`grendel.models.hnl.production`); the
lifetimes and visible final states come from the decay templates, not
from HNLCalc.

### FairShip — the SHiP Collaboration — LGPL-3.0-or-later

The SHiP Collaboration, *FairShip*, https://doi.org/10.5281/zenodo.18020628
(concept DOI, all versions); `ShipSoft/FairShip`; Copyright CERN for the
benefit of the SHiP Collaboration. Cite the DOI wherever an HNL width,
lifetime or decay table from FairShip enters a result. Used only for the v1
HNL templates and the width-uncertainty band; no v2 figure needs it.

Only its `python/` directory is fetched, at the revision in `PINS.json`
(4 March 2026), which is on `main` between releases 25.12 and 26.03. Of the
files used, `hnl.py`, `shipunit.py` and `alpha_s.dat` are byte-identical to
release 26.03, and the differences in `readDecayTable.py` and
`pythia8_conf_utils.py` are lint-only, so release 26.03 is the archived
version this pin corresponds to. The modules supply the flavour-dependent
HNL widths and the Pythia 8 decay-channel configuration behind the FairShip
decay templates of the first published curve and the width-band table
(`grendel.models.hnl.templates.fairship_decay`, `grendel.band.hnl.width_band`).
They are used unmodified: the pipeline puts the directory on `sys.path` and
imports ROOT itself before loading them. Needs PyROOT with Pythia 8
(`ROOT.TPythia8`); FairShip itself is not built.

**Code in this repository derived from FairShip.** Two files are adaptations
of FairShip sources rather than callers of them, and carry FairShip's licence
(SPDX headers in each; text in `LICENSES/LGPL-3.0-or-later.txt`):

| File | Derived from | Licence |
|---|---|---|
| `grendel/models/hnl/templates/fairship_decay.py` | `python/readDecayTable.py::addHNLdecayChannels` (in `_add_hnl_decay_channels`) and the HNL particle and decay-mode command block of `python/pythia8_conf.py::add_hnl` (in `FairShipDecayBackend._configure_pythia`) | LGPL-3.0-or-later |
| `grendel/models/hnl/data/fairship_decay_selection.conf` | `python/DecaySelection.conf`, every channel switched on | LGPL-3.0-or-later |

`_add_hnl_decay_channels` exists because `addHNLdecayChannels` writes through
a `ROOT.HNLPythia8Generator`, which needs a FairShip build; this analysis
drives a bare `ROOT.TPythia8` instead. The physics — the channel list, the
`BR/2` Majorana split, the charge-conjugate modes — is FairShip's.

Everything else that touches FairShip only imports it at run time and stays
BSD-3-Clause, which LGPL-3.0-or-later permits: `width_band.py` calls
`hnl.HNLbranchings` and applies its own floor/cap/envelope rule to the
resulting widths, and the ALP template generator borrows only this
repository's own ROOT plumbing. The tables produced by running FairShip
(`width_band_delta.csv`, the cached decay templates) are output, not
derived source.

### HeavyN UFO model `SM_HeavyN_CKM_AllMasses_LO` — R. Ruiz and collaborators — FeynRules model database

Version 2.2 (R. Ruiz, 27 June 2016), downloaded from the FeynRules HeavyN
model page, archive SHA-256 in `PINS.json`. D. Alva, T. Han, R. Ruiz,
arXiv:1411.7305; C. Degrande, O. Mattelaer,
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

### Pythia 8.317 — C. Bierlich et al. — GPL-2.0-or-later

`Pythia8/releases` tag `pythia8317`. Bierlich et al., SciPost Phys.
Codebases 8 (2022), arXiv:2203.11601. The only Pythia needed: exHad requires
it, it makes the tracked kaon spectrum
(`grendel.models.hnl.production.make_kaon_spectrum`) and it showers the BC5
Higgs events (`grendel.models.scalar.higgs_shower`, the committed (pT, y)
histogram). Both are C++ drivers, so no Pythia Python bindings are needed. Up
to v2.2 the kaon spectrum was made with 8.315 (0.5 % in <n_K>) and the Higgs
histogram with 8.312 (mean pT and |y| within 0.2 % of 8.317 in every mode);
the v1 FairShip templates also used 8.312.

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
  O. Mikulenko, K. Bondarenko, arXiv:2305.13383; BSD-3-Clause) v1.4, the revised phenomenology of arXiv:2501.04525; its
  `PROVENANCE.md` records the export.
- `grendel/models/alp/data/alpinist/` — digitised GKOZ (arXiv:2310.03524)
  width curves from ALPINIST (J. Jerhot, B. Döbrich, F. Ertas, F. Kahlhoefer,
  T. Spadaro, arXiv:2201.05170; BSD-3-Clause; licence in
  `curves/sources/alpinist/LICENSE`), retained for audit; its `COMMIT_SHA.txt`
  also pins the ALPINIST source of `tools/compute_cbs_alpinist.py`.
- `grendel/models/scalar/data/winkler_widths.csv` — dark-scalar hadronic
  widths digitised from M. W. Winkler, Phys. Rev. D 99 (2019) 015018,
  arXiv:1809.01876, Fig. 4.
- `grendel/models/scalar/data/eventcalc_1809/` — the dark-scalar lifetime
  and branching ratios of the same model to 63 GeV, the EventCalc tables
  shipped with exHad (Kryshtal, Ovchynnikov, arXiv:2609.16104, at the pinned
  commit; EventCalc / SensCalc, Ovchynnikov et al., arXiv:2305.13383;
  BSD-3-Clause, notice in the directory). Used by BC5 above 7.5 GeV; its
  `PROVENANCE.md` records the mapping.
- `grendel/models/scalar/data/senscalc_quartic/` — BR(B+ → X_s S S) and
  BR(B_s → S S) of the quartic coupling from SensCalc v1.4 (BSD-3-Clause,
  `LICENSE` in the directory), the phenomenology of Boiarska et al.,
  arXiv:1904.10447; `PROVENANCE.md` records the normalisation.
- `grendel/models/scalar/data/higgs_pt_y_14TeV.csv.gz` — Higgs (pT, y)
  histograms per production mode generated here with MadGraph5_aMC@NLO
  v3.6.6 (LO, NNPDF4.0 NLO) and Pythia 8 (its header records the versions
  and seeds), normalised in the code to the LHC Higgs Cross Section Working
  Group YR4 14 TeV cross sections (arXiv:1610.07922).
- `curves/data/` — reference curves digitised from the figures of the
  papers named in each file's header (PBC reports arXiv:1901.09966,
  arXiv:2310.17726, arXiv:2505.00947; SHiP, FASER2, CODEX-b, ANUBIS,
  MATHUSLA40, DarkQuest, GKOZ, FIPs 2022, the BBN bound of Chen and Zhang,
  arXiv:2410.07343), the HNLimits workbook, and
  the ALPINIST NA62 bound; see `curves/README.md` and `curves/PROVENANCE.md`.
- `grendel/production/data/fonll/central/` — meson grids produced with the
  software above (`grendel/production/data/fonll/PROVENANCE.md`); cite
  FONLL and NNPDF4.0 (arXiv:2109.02653) when using them.
- `grendel/geometry/` — the detector geometry and reconstruction shared
  with the Higgs-portal study of the same paper; `grendel/geometry/PROVENANCE.md`.
  Since October 2026 the tunnel is UP542 as built: `grendel/geometry/data/up542/`
  is an unchanged mirror of `results/geometry/` in
  https://gitlab.cern.ch/grendel/grendel-geometry (CERN-internal).
