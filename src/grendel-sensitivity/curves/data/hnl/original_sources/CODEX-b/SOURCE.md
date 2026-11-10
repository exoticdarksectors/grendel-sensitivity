# CODEX-b — original-source provenance (way 4)

**Citations:** arXiv:1911.00481, arXiv:2203.07316

**Figures containing the contours:** `figures/HNL_{Ele,Mu,Tau}.pdf` of
1911.00481v2 (`fig:{Ne,mu,tau}_combined`) — the "Combined" single-flavour
panels.

**Note:** No public machine-readable data, but **way-4 is done** (2026-08-28).
`curves.tools.vector_codexb_hnl` vector-extracts the CODEX-b contour from the
EoI arXiv source's own three per-flavour panels — all of Ue/Umu/Utau, unlike
the muon-only PBC BC7 redraw — and converts the EoI's Dirac convention to
Majorana. The 2025 CODEX-b ESPP update (arXiv:2505.05952) redraws the EoI's
Dirac contour unchanged in its Fig. 9 (muon and tau panels, "modified from" the
EoI; checked 2026-10-25), so 1911.00481 remains the citable source.

**Data files (way 4):**
- `CODEX-b_2019_Ue.dat`   `96396c25a715f39aedad5c88dbbb3337dd6998379600051d2ab60611a2cb2a9c cef2ed82543c3642b24b70108fd0944a88f94f223b4c271df77e17e81a323f16 f503b91e7918884d20e700b48f354690efb6d2e7ea548ecf3ac7c392d4059719 `
- `CODEX-b_2019_Umu.dat`  ``
- `CODEX-b_2019_Utau.dat` ``

Columns are `mass_GeV  U2_majorana  U2_dirac_raw`; NaN marks a vertex one of
the two contours does not have (drop a column's NaN rows to get its polyline).

**Source figure PDFs:** `sources/arXiv-1911.00481/HNL_{Ele,Mu,Tau}.pdf`, copied
from the arXiv v2 source bundle's `figures/` directory.

**Method:** CODEX-b is isolated by stroke colour `RGB(1.00,0.55,0.00)`,
confirmed against the legend swatch in each panel rather than assumed. Axes are
calibrated from **tick geometry, never from tick-label text bounding boxes** —
a superscripted `10^-6` has its mantissa glyph centred below the tick, which is
exactly the defect that biased the BC4 CODEX-b extraction by 8.6% (see
`PROVENANCE.md`). All three panels share x 150.3947 pt/decade,
y 32.1590 pt/decade, frame 0.2–50 GeV and 10^-2..10^-10.

**Nature:** the EoI caption states **Dirac** HNLs; GRENDEL and the HNLimits
bounds are Majorana. Column 2 follows the CODEX-b authors' prescription
(S. Knapen, e-mail of 2026-10-05, since 2026-10-06): lower-edge vertices
`/ sqrt(2)`, lifetime-limited upper-edge vertices `/ 2`, the tip at the
Dirac tip mass (Ue 3.10, Umu 3.00, Utau 3.20 GeV). With equal production and
`Gamma_M = 2 Gamma_D` at fixed `U^2`, `N_M(U^2) = N_D(2 U^2)/2`: on the lower
edge `N ~ U^4` and the shift is exactly `1/sqrt(2)`; on the upper edge `U^2/2`
keeps the lifetime but drops the halved production, so it is the outer bound
of the Majorana island (the exact edge lies lower by a relative
`~ln2/(a U^2)`, a few per cent). The frame vertex at `U^2 = 1e-2` has no
Majorana value, and the curve stays open at the top. Column 3 is the
unconverted value as drawn. Before: per-column `N_D(x) = C x^2 exp(-a x)`
model (2026-09-28 to 2026-10-06), and HNLimits' uniform
`dirac_to_majorana_dic["BD"] = 1/sqrt(2)` on every vertex (v1).

**Sanity check:** best reach as drawn is `1.88e-08` / `1.99e-08` / `8.38e-07`
for Ue/Umu/Utau, matching the EoI's own stated ~1e-8 for e, µ and ~1e-6 for τ;
after conversion `1.33e-08` / `1.40e-08` / `5.93e-07` (at 1.5 / 1.5 / 1.2 GeV).

**arXiv PDFs cached at:**
- 1911.00481: cached -> <work dir>/curves/original_sources_pdfs/1911.00481.pdf
- 2203.07316: cached -> <work dir>/curves/original_sources_pdfs/2203.07316.pdf

**Note on regeneration:** this file is hand-maintained. The script
`curves.tools.original_sources` marks CODEX-b (and ANUBIS) `hand_maintained`
and leaves their `SOURCE.md` alone; `curves.tools.vector_codexb_hnl` rewrites
the data files, whose hashes above must then be updated by hand.
