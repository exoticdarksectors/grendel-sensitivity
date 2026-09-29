# Third-party software and data

External programs are fetched by `third_party/fetch.py` at the revisions in
`third_party/PINS.json` and keep their own licences.

| Package | Licence | Used for | Changes |
|---|---|---|---|
| exHad — Kryshtal, Ovchynnikov, arXiv:2609.16104 | BSD-3-Clause | exclusive hadronic decay templates (HNL, BC4, BC10) | none |
| HNLCalc — Feng, Hewitt, Kling, La Rocco, arXiv:2405.07330 | none published | HNL production rates | fork `leoredi/HNLCalc`: an unused import dropped; a dead duplicate D_s → K⁰ form factor removed; the τ → P N Källén factor and the b-baryon three-body mixing power corrected |
| FairShip — SHiP Collaboration, doi:10.5281/zenodo.18020628 | LGPL-3.0-or-later | HNL widths and decay tables | none |
| FONLL 1.3.3 — Cacciari, Frixione, Nason; cite JHEP 05 (1998) 007 and JHEP 03 (2001) 006 | its authors' terms | heavy-flavour grids | `grendel/production/fonll_grids/patches` |
| MadGraph5_aMC@NLO 3.6.6 — arXiv:1405.0301 | MadGraph licence | W/Z → ℓN and prompt-τ samples | none |
| HeavyN UFO `SM_HeavyN_CKM_AllMasses_LO` — arXiv:1411.7305, 1602.06957, 1812.08750 | FeynRules model database | the MadGraph HNL model | none |
| Pythia 8.315 and 8.317 — arXiv:2203.11601 | GPL-2.0-or-later | the charged-kaon spectrum (8.315); exHad (8.317) | none |

`grendel/models/hnl/templates/fairship_decay.py` and
`grendel/models/hnl/data/fairship_decay_selection.conf` follow FairShip
sources and are LGPL-3.0-or-later (SPDX headers; licence text in
`LICENSES/`). The rest of the repository is BSD-3-Clause.

## Data

- `grendel/production/data/fonll/central/` — FONLL grids with NNPDF4.0
  (arXiv:2109.02653).
- `grendel/models/alp/data/senscalc_2501/` — ALP widths and branching ratios
  from SensCalc (arXiv:2305.13383, BSD-3-Clause) for the model of
  arXiv:2501.04525.
- `grendel/models/scalar/data/winkler_widths.csv` — dark-scalar widths from
  M. W. Winkler, arXiv:1809.01876, Fig. 4.
- `grendel/models/scalar/data/senscalc_quartic/` — B → K S S, b → X_s S S and
  B_s → S S rates from SensCalc (BSD-3-Clause, `LICENSE`), in the conventions
  of Boiarska et al., arXiv:1904.10447.
- `grendel/models/scalar/data/eventcalc_1809/` — dark-scalar branching ratios
  and lifetimes above 7.5 GeV from EventCalc (BSD-3-Clause,
  `LICENSE-EventCalc`).
- `grendel/models/scalar/data/higgs_pt_y_14TeV.csv.gz` — Higgs (pT, y) at
  14 TeV from MadGraph and Pythia 8 over all production modes, normalised to
  the LHC Higgs Cross Section Working Group, arXiv:1610.07922.
- `grendel/models/hnl/production/data/kaon_softqcd_spectrum.npz` — Pythia
  8.315 SoftQCD charged-kaon spectrum.
- The BC10 b → s coupling and B → K⁽ⁱ⁾ form factors follow GKOZ
  (arXiv:2310.03524) as implemented in ALPINIST (arXiv:2201.05170).
- `grendel/geometry/grendel_geometry.py` and `reco_common.py` come from
  `exoticdarksectors/llpatcolliders`.
