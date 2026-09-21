# CODEX-b — original-source provenance (way 4)

**Citations:** arXiv:1911.00481, arXiv:2203.07316

**Figures containing the contours:** `figures/HNL_{Ele,Mu,Tau}.pdf` of
1911.00481v2 (`fig:{Ne,mu,tau}_combined`) — the "Combined" single-flavour
panels.

**Note:** No public machine-readable data, but **way-4 is done** (2026-08-28).
`curves.tools.vector_codexb_hnl` vector-extracts the CODEX-b contour from the
EoI arXiv source's own three per-flavour panels — all of Ue/Umu/Utau, unlike
the muon-only PBC BC7 redraw — and converts the EoI's Dirac convention to
Majorana. The 2025 CODEX-b ESPP update (arXiv:2505.05952) carries no HNL
contour of its own and defers to this EoI, so 1911.00481 remains the citable
source.

**Data files (way 4):**
- `CODEX-b_2019_Ue.dat`   `df35b8c13c25ff132442035c8389425a996120dffe4786fc55769f85497c8b05`
- `CODEX-b_2019_Umu.dat`  `c29c152a2ac7a7816b3276f28465be3d5d1ea1a0fed9cb76dd2ad04c9356c393`
- `CODEX-b_2019_Utau.dat` `46cf181ccb9fc7235d37b84ec2bfec2bb9c3e9840edcde7c2dbf2678ca489628`

Columns are `mass_GeV  U2_majorana  U2_dirac_raw`.

**Source figure PDFs:** `sources/arXiv-1911.00481/HNL_{Ele,Mu,Tau}.pdf`, copied
from the arXiv v2 source bundle's `figures/` directory.

**Method:** CODEX-b is isolated by stroke colour `RGB(1.00,0.55,0.00)`,
confirmed against the legend swatch in each panel rather than assumed. Axes are
calibrated from **tick geometry, never from tick-label text bounding boxes** —
a superscripted `10^-6` has its mantissa glyph centred below the tick, which is
exactly the defect that biased the BC4 CODEX-b extraction by 8.6% (see
`PROVENANCE.md`). All three panels share x 150.3947 pt/decade,
y 32.1590 pt/decade, frame 0.2–50 GeV and 10^-2..10^-10.

**Nature:** the EoI caption states **Dirac** HNLs. Column 2 applies HNLimits'
factor for a production × decay limit, `dirac_to_majorana_dic["BD"] =
1/sqrt(2)`, matching the Majorana GRENDEL curve and HNLimits bounds it is drawn
against. Column 3 is the unconverted value as drawn.

**Sanity check:** best reach as drawn is `1.88e-08` / `1.99e-08` / `8.38e-07`
for Ue/Umu/Utau, matching the EoI's own stated ~1e-8 for e, µ and ~1e-6 for τ.

**arXiv PDFs cached at:**
- 1911.00481: cached -> <work dir>/curves/original_sources_pdfs/1911.00481.pdf
- 2203.07316: cached -> <work dir>/curves/original_sources_pdfs/2203.07316.pdf

**Note on regeneration:** this file is normally written by
`curves.tools.original_sources`. It was hand-updated on 2026-08-28 to record the
completed way-4 extraction; the matching `note` text in that script was updated
at the same time, so re-running it reproduces the same claims.
