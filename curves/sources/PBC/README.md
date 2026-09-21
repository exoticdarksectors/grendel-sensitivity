# PBC report source bundles

Unpacked arXiv source bundles of the three Physics Beyond Colliders (PBC)
documents that define and track the BC benchmark models. The bundles
themselves are gitignored (re-downloadable from arXiv); this README records
what each one is and which benchmarks it covers, because the split between
them is a recurring source of confusion (reviewed against the sources
2026-07-22).

## The three documents

- `arXiv-1901.09966v2/` — Beacham et al., *Physics Beyond Colliders at CERN:
  BSM Working Group Report* (2019), J. Phys. G 47 (2020) 010501.
  **This is the document that defines the benchmarks BC1–BC11** (Section 2.1
  of the published layout; `PBC_BSM_document_v1.tex` here). Everything else
  cites back to it (`Beacham:2019nyx`).
- `arXiv-2310.17726v1/` — Ahdida et al., *Post-LS3 Experimental Options in
  ECN3*, CERN-PBC-REPORT-2023-003 (2023). Sensitivity update for the
  beam-dump-oriented benchmarks: **BC1, BC2, BC4, BC6, BC8, BC9, BC10**
  (`9-PhysicsPotential.tex:172` and Figures BC1–BC10). The BC10 figure is
  restricted to `m_a < m_eta'` "to avoid large theoretical uncertainties at
  higher masses".
- `arXiv-2505.00947v2/` — Antel et al., PBC strategy update (2025).
  Sensitivity update for the complementary, LHC-adjacent subset:
  **BC3, BC5, BC7** plus generic `h -> XX` / `h -> A'A'` showcase plots.
  It covers only benchmarks *not* detailed in the ECN3 report — stated
  explicitly at `5-PhysicsReach.tex:339`: "The former BCs have been
  presented in detail in [Ahdida:2867743]; therefore we focus here on the
  latter."

Together the two update reports cover BC1–BC10; **BC11** (gluon-coupled ALP)
appears in neither update, only in the 2019 definitions.

## Benchmark coverage map

| BC   | Model                          | 2019 defines | ECN3 2023 fig. | 2025 update fig. |
|------|--------------------------------|:---:|:---:|:---:|
| BC1  | Dark photon, visible           | yes | yes | —   |
| BC2  | Dark photon -> invisible DM    | yes | yes | —   |
| BC3  | Millicharged particles         | yes | —   | yes |
| BC4  | Dark scalar, theta only        | yes | yes | —   |
| BC5  | Dark scalar + h->SS production | yes | —   | yes |
| BC6  | HNL, e dominance               | yes | yes | —   |
| BC7  | HNL, mu dominance              | yes | —   | yes |
| BC8  | HNL, tau dominance             | yes | yes | —   |
| BC9  | ALP–photon                     | yes | yes | —   |
| BC10 | ALP–fermion                    | yes | yes | —   |
| BC11 | ALP–gluon                      | yes | —   | —   |

Practical consequences for the digitization scripts here:

- The dark-scalar "currently excluded" region and competitor projections are
  taken from the 2025 report's **BC5** figure (`pbc_bc5_gray_final.pdf`) —
  the mixing-angle limits shown there apply to BC4 unchanged.
- The 2025 report has **no BC10 figure**; the BC10 excluded region is
  extracted from GKOZ (arXiv:2310.03524) instead. A PBC BC10 sensitivity
  figure *does* exist — in the ECN3 report
  (`Figures/9-PhysicsPotential/BC10.pdf`) — it is simply not in the 2025
  update.

## Single-coupling structure (verbatim, 2019 report)

All GRENDEL single-portal benchmarks have production and lifetime locked
together by one parameter, by construction of the benchmark:

- BC4: "in this model we assume lambda = 0, and all production and decay are
  controlled by the same parameter theta" (`PBC_BSM_document_v1.tex:623`).
- BC6–BC8: "all production and decay can be determined as function of
  parameter space (m_N, |U_l|^2)".
- BC10: "all phenomenology (production and decay) can be determined as
  functions on {m_a, f_l^-1, f_q^-1}. Furthermore, for the sake of
  simplicity, we take f_q = f_l" (`PBC_BSM_document_v1.tex:710`).

By contrast, the `h -> SS` presentation (BC5-adjacent) deliberately treats
BR(h->SS) and c*tau as independent free parameters at fixed mass — the 2025
report's showcase figures (`fig:pbc_h_to_scalar`, `fig:pbc_h_to_dark_photon`)
use this convention, which is why exotic-Higgs sensitivities are quoted at a
fixed representative mass (e.g. 15 GeV) rather than as closed islands in a
(mass, coupling) plane.
