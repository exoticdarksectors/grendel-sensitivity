# Source figures of the reference curves

The extraction tools in `curves/tools` read the published figures listed
here. None is tracked; each is an arXiv source-bundle file, downloaded from
`https://arxiv.org/e-print/<id>` and unpacked so that the paths below exist
under this directory (or under `GRENDEL_CURVES_SOURCES`). The tracked
outputs in `curves/data` are current, so nothing needs to be downloaded to
render the figures; the sources are only needed to *re-extract* a curve.

| Path under `curves/sources` | Paper | Used by |
|---|---|---|
| `PBC/arXiv-2505.00947v2/Figures/5-PhysicsReach/BSM_benchmarks/pbc_bc7_gray_final.pdf` | Antel et al., PBC strategy update 2025 | `vector_bc7`, `raster_bc7` |
| `PBC/arXiv-2505.00947v2/Figures/5-PhysicsReach/BSM_benchmarks/pbc_bc5_gray_final.pdf` | idem | `vector_bc4` (existing exclusion; FASER2/SHiP BC5 references) |
| `arXiv-2504.06692/BC4.pdf` | SHiP Collaboration 2025 (panel present in the arXiv source only) | `vector_bc4` (SHiP) |
| `arXiv-1911.00481/B_decay_benchmark.pdf` | CODEX-b Expression of Interest | `vector_bc4` (CODEX-b) |
| `arXiv-1911.00481/HNL_{Ele,Mu,Tau}.pdf` | idem | `vector_codexb_hnl` |
| `arXiv-2310.03524v3/LHCb-constraints.pdf` | GKOZ, fermion-coupled ALP phenomenology | `vector_bc10` (LHCb recast) |
| `arXiv-2501.04525/parameter-space-ALP-fermion.pdf` | Ovchynnikov–Zaporozhchenko, PRD 112 015001 | `vector_bc10` (SHiP, DarkQuest, past beam dumps) |
| `arXiv-2305.01715v1/Light_DM/summary_plots_LDM/ALPs_fermions_fips2022_text_logo.pdf` | FIPs 2022 workshop report | `vector_bc10_fips2022` |
| `arXiv-2606.26862v1/img/HNL_{electron,muon}_ANUBIS_Combination_total_sigEff0p5.pdf` | ANUBIS HNL projection 2026 | `data/hnl/original_sources/ANUBIS` (extracted with the `pdfio` helpers; see its `SOURCE.md`) |
| `arXiv-2512.13011v1/8_sensitivity_VvSq_vs_mass_DandB_merged_LHCandSHiP.pdf` | Wang–Zhang, ANUBIS meson production | idem |
| `alpinist/` (tracked) | ALPINIST `Figures/Bound_data/gY`, commit in `COMMIT_SHA.txt`, BSD-3 | `vector_bc10` (NA62 band) |

The HNLimits curves come from the installed `HNLimits` package
(`pip install -e ".[curves]"`), not from here; the workbook that carries
their metadata is tracked as `curves/data/hnl/local_HNL_database.xlsx`.

`PBC/README.md` records which of the three PBC reports defines and which
one updates each benchmark, since the split is a recurring source of
confusion.
