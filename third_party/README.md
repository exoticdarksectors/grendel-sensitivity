# External software

The pipeline drives several programs and models that are not part of this
repository. `fetch.py` checks each one out here, from its own upstream, at
the revision recorded in `PINS.json`:

```
python third_party/fetch.py --list          # what, from where, at which pin, and what to do next
python third_party/fetch.py exHad HNLCalc   # fetch some
python third_party/fetch.py --all
python third_party/fetch.py --verify        # is everything present at its pin?
```

`GRENDEL_THIRD_PARTY_DIR` relocates the checkouts; the package resolves
them through `grendel.io.paths.third_party_dir()`. The checkouts themselves
are ignored by git; only `fetch.py`, `PINS.json` and this file are tracked.

| Name | What it is | Used for | After fetching |
|---|---|---|---|
| `exHad` | hadronic decays of GeV-scale feebly interacting particles (Kryshtal, Ovchynnikov) | the HNL, BC4 and BC10 rest-frame decay templates | `python tools/setup.py --pythia8-dir <pythia8317>` in its own venv (`GRENDEL_EXHAD_PYTHON`) |
| `HNLCalc` | HNL production and decay rates (Feng, Hewitt, Kling, La Rocco) | meson, baryon and tau production branching ratios | none (pure Python) |
| `FairShip` | the SHiP software's HNL modules (`python/` only) | the FairShip decay templates of the first published curve; the width-band table | a PyROOT + Pythia8 environment |
| `SM_HeavyN_CKM_AllMasses_LO` | the HeavyN FeynRules/UFO model (Ruiz et al.) | W/Z -> l N production with MadGraph | none; MadGraph converts it on first import |
| `fonll` | FONLL 1.3.3 heavy-quark production | regenerating the meson grids | `python -m grendel.production.fonll_grids.install --build` |
| `MG5_aMC_v3_6_6` | MadGraph5_aMC@NLO 3.6.6 | W/Z and tau production of HNLs | none beyond its own first-run setup; needs LHAPDF (`GRENDEL_LHAPDF_CONFIG`) |
| `pythia8315` | Pythia 8.315 | the K -> l N kaon spectrum (reproduces the tracked spectrum only at 8.315) | `./configure && make` |
| `pythia8317` | Pythia 8.317 | exHad's event generation (exHad requires exactly 8.317) | `./configure && make`, then exHad's `tools/setup.py` |

Licences, citations and the modifications carried by the pinned revisions
are in `THIRD_PARTY.md` at the repository root.
