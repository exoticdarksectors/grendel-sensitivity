# Third-party code, software and data

Two pieces of code are in this repository. Everything else the pipelines use is downloaded at a fixed version, by
`env/bootstrap.sh` or by the first step of a pipeline. It keeps its own licence and is not redistributed here.

## In this repository

| Path | What | Licence |
|---|---|---|
| `src/grendel-sensitivity/` | the analysis code grendel-sensitivity v2.5.0: `git archive` of commit 8802f50 of its development repository, unmodified | BSD-3-Clause and LGPL-3.0-or-later: its `LICENSE`, `LICENSES/` and `THIRD_PARTY.md` (data derived from published work, with their citations) |
| `src/llpatcolliders-hnl-06e5a18/` | `hnl/` of [exoticdarksectors/llpatcolliders](https://github.com/exoticdarksectors/llpatcolliders) at 06e5a18 (2026-06-23; public through its pull request #14, branch `hnl-production-mesons-ew`), the June 2026 HNL production code that some HNL channels still use; unmodified apart from the added `COMMIT` and `LICENSE` | BSD-3-Clause, (c) 2026 Federico Leo Redi (`LICENSE` in that folder) |
| `reference/paper_figures/` | the four PNG figures of Sec. IV of arXiv:2609.00152v2 (Figs. 4-7), as in the paper | (c) the authors of the paper |
| everything else | the reproduction: environment, pipelines, figures, results | BSD-3-Clause, (c) 2026 Federico Leo Redi (`LICENSE`) |

## Downloaded at fixed versions

| Software or data | Version | Licence | Fetched by |
|---|---|---|---|
| exHad (V. Kryshtal, M. Ovchynnikov) | e41b07b | BSD-3-Clause | `env/bootstrap.sh` (`third_party/fetch.py`, pins in `src/grendel-sensitivity/third_party/PINS.json`) |
| HNLCalc (fork with six fixes, upstream laroccod/HNLCalc) | leoredi/HNLCalc b1d36c6 | none published upstream | same |
| HeavyN UFO model `SM_HeavyN_CKM_AllMasses_LO` (FeynRules database) | sha256 05dab629... | FeynRules model database | same |
| FONLL 1.3.3 (alisw mirror) | c7086e4 | distributed by its authors | same |
| MadGraph5_aMC@NLO | 3.6.6 (fb1c2f2) | UoI-NCSA | same |
| Pythia | 8.317 (83a5c99) | GPL-2.0-or-later | same |
| conda-forge toolchain (Python 3.11, gcc/gfortran, cmake, LHAPDF 6.5, ...) | exact builds in `env/conda-linux-64.lock` | each package's own | `env/bootstrap.sh` |
| Python packages (numpy, pandas, matplotlib, HNLimits 1.2.0 (MIT), particle, ...) | exact versions in `env/requirements-venv.lock.txt`, `env/requirements-exhad.lock.txt` | each package's own | `env/bootstrap.sh` |
| `particle` 0.26.2 (BC4 production only, as published) | 0.26.2 | BSD-3-Clause | `pipelines/scalar/00_setup.sh` |
| PDF set NNPDF40_nlo_as_01180 (LHAPDF) | DataVersion 1, tarball sha256 a2ef4869... | NNPDF | `env/bootstrap.sh` |
| PDF set cteq66 (LHAPDF) | files checked by sha256 | CTEQ | `pipelines/fonll/00_fetch_pdfs.sh` |
| Wolfram Engine | 15.0.0 | proprietary; the free Wolfram Engine for Developers licence. Each user activates it with their own Wolfram ID; it cannot be redistributed | `pipelines/senscalc/00_install_wolfram.sh` |
| SensCalc (M. Ovchynnikov et al.) | eb2c456 (and v1.4, 4f7210b) | BSD-3-Clause | `pipelines/senscalc/01_fetch_sources.sh` |
| FeynCalc | 10.2.1 | GPL-3.0 | same |

`src/grendel-sensitivity/CITATION.cff` and its `THIRD_PARTY.md` give the paper to cite for each tool.
