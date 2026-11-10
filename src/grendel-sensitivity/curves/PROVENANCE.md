# Figure provenance and conventions

The three sensitivity figures of arXiv:2609.00152 are rendered by
`python -m curves.plot paper` from a results tree holding the central
`N_sig >= 3` scan and the same scan re-solved at `N_sig >= 10`; the
published (v1) tables are tracked as that tree in `data/grendel_v1`
(`python -m curves.plot paper --grendel-dir curves/data/grendel_v1 --benchmark hnl bc4 bc10`):

| Paper figure | File | Layers |
|---|---|---|
| Fig. 4 | `bc4_grendel_paper.png` | current exclusion (PBC 2025 BC5 figure), SHiP and FASER2 3 ab⁻¹ (ESPP 2026 Briefing Book), CODEX-b 300 fb⁻¹, GRENDEL |
| Fig. 5 | `hnlimits_grendel_paper.png` | HNLimits experimental bounds, single-Majorana BBN bound (Chen & Zhang), SHiP, FASER2, CODEX-b, GRENDEL; one panel per flavour |
| Fig. 6 | `bc10_grendel_paper.png` | LHCb `B → K μμ` recast (GKOZ), past beam dumps, SHiP, GRENDEL |
| — | `bc5_grendel_paper.png` | not in the manuscript: current exclusion, CODEX-b, MATHUSLA40, SHiP (PBC 2025 BC5 figure), ANUBIS (arXiv:2510.26932, conservative background), GRENDEL |

The PNGs published in arXiv:2609.00152v1 have SHA-256

```
bc4_grendel_paper.png       7da7e47e5c741c0580f398f00d4dcfbc9cb058928c29e8bc0ea07424e7e9ec21
hnlimits_grendel_paper.png  dcad42ae1fc6aa206a289b8f27bae6368b8fb86af94aa593fdc6d41f6d71e71a
bc10_grendel_paper.png      d73ee0f8423e0a6937600c203d6e6b66c0ac899520beced66bc750fafcce68cd
```

and the renderer up to 2.2.1, with the v1 comparison layers, reproduces
them byte for byte from `data/grendel_v1` on macOS with matplotlib 3.10.9 and HNLimits 1.2.0
(Python 3.11 / numpy 2.4.6 at the v1 freeze; re-checked with Python 3.13 /
numpy 2.5.3). A Linux build of the same versions (the `python:3.11`
container) draws identical pixels but compresses the PNG differently, so
compare decoded pixels across platforms. Four comparison layers have since
been corrected and, in 2.2.2, every panel restyled (one colour and dash
pattern per experiment, GRENDEL solid, thinner lines, one legend order), so
every panel differs from v1 while the GRENDEL curves do not; the current
renderer gives, on macOS,

```
bc4_grendel_paper.png       583ac4c5626e4445169c5fc27cf0acc47fac83ef20e039ce37e1ac74454f5b86
hnlimits_grendel_paper.png  284631ff2f02a4ca1aa4b9b71c8068db5916ef2b4fbdae0f2d2f16897a92a2d0
bc10_grendel_paper.png      fda72747d9f6ed13943df7ae82281a613cfa43f6198fb96e867b315944a3b464
```

| Figure | Layer changed relative to v1 |
|---|---|
| Fig. 4 | SHiP and FASER2 from the BC4 panel of the ESPP 2026 Briefing Book (v1 drew the 2504.06692 SHiP trace and no FASER2); luminosities in the legend; frame floor 1e-14, the Briefing Book's |
| Fig. 5 | cosmology: the single-Majorana BBN bound of Chen & Zhang replaces Sabti et al. × 0.41, legend "BBN bound (single Majorana HNL)"; CODEX-b: per-column Dirac → Majorana conversion replaces a uniform 1/√2 |
| Fig. 6 | past beam dumps clipped to their source frame, 0.6–2.0 GeV (v1 filled them from 0.555 GeV) |

The PNG is the deterministic artefact; matplotlib stamps a creation date
into the PDF written alongside it, so PDF hashes change on every run.

The paper profiles select layers and the final-size layout only. They
call the same loaders, metadata transformations, clipping rules and
central-curve renderer as the talk layout, so the two cannot diverge; the
diagnostic layout reads the raw HNLimits copies with its own renderer.

## Two signal thresholds

Every panel draws the `N_sig >= 3` island as the filled red band and the
`N_sig >= 10` island as a solid dark-red contour nested inside it: GRENDEL
alone is drawn solid, and every projection is broken, one colour and dash
pattern per experiment on every panel (`common.PROJECTION_STYLES`). The panels and
the paper write `N_sig > 3`; the solver compares `N >= N_THRESHOLD`
(`grendel.constants`). `N` is a continuous expected yield, so both
spellings select the same contour.

The two are not independent calculations. The geometry ray-cast, the
acceptance Monte Carlo and the coupling scan are threshold independent;
only the island solver consumes the threshold, and `run_scan --thresholds
3 10` extracts both islands from one yield scan per mass. On the
production-limited lower edge both the production rate and the decay
probability scale with the coupling, so the `N>=10 / N>=3` edge ratio has
an analytic limit: BC4 (`N ~ sin⁴θ`) predicts `(10/3)^(1/2) = 1.826`, BC10
(`N ~ (1/f)⁴`) predicts `(10/3)^(1/4) = 1.351`; the published ratios bottom
out at 1.826 and 1.353, at the lowest masses.

Island tips are closed by interpolating the yield-threshold crossing in
log yield and log coupling between adjacent grid points
(`insert_threshold_tips`, `append_high_mass_tip`); a curve solved at a
different threshold closes at that threshold.

## HNL panels (BC6–BC8)

The experimental bounds are HNLimits' curves transformed with the unit,
confidence-level and Dirac/Majorana metadata of its workbook
(`curves/data/hnl/local_HNL_database.xlsx`, SHA-256
`aa9346350fd2354af29b71c3a4c6d2da1c5454abd6cef3b780dd1dd858f56783`); the
raw `.dat` files alone lose those transformations. The drawn regions,
cosmology included, are written out by `curves.tools.hnlimits_process` and
checked in CI against the tracked copies.

Three local decisions depart from what HNLimits would draw:

- **Cosmology is a single-Majorana BBN bound.** HNLimits' `cosmo` row is
  the BBN bound of Sabti et al. (arXiv:2006.07387), simulated for a Dirac
  HNL, i.e. two degenerate Majorana HNLs, of which their Sec. 5.2 says "the
  bounds do not apply to models with single Majorana HNLs". The paper's
  benchmark is a single Majorana HNL, so the layer is replaced (v1 drew it,
  per flavour, times HNLimits' Dirac and 2σ factors 1/2 × √(2.71/4) = 0.41,
  a second conversion of a curve already per Majorana state). Boyarsky et
  al. (arXiv:2008.00749) also assume a quasi-Dirac pair; HNLimits ships no
  other BBN curve, the PBC 2025 BC7 figure draws none and the ECN3 BC6/BC8
  figures an uncited one. Chen and Zhang (arXiv:2410.07343, PRD 111 (2025)
  123024) compute exactly the benchmark: one Majorana HNL mixing with one
  flavour, abundance and decays both fixed by `(m_N, |U|^2)`, excluded at
  95% CL. `COSMO_PER_FLAVOUR_FILES` draws their per-flavour region as
  published, no CL or Dirac/Majorana rescaling
  (`data/hnl/cosmology/`, vector-extracted by `curves.tools.vector_bbn_hnl`
  from their Figs. 1 and 8). Their calculation stops at 1.2 GeV, leaving
  1–10 GeV as a gap they expect to interpolate between the computed
  regions; the upper edge is below the `1e-10` floor by 1.01 GeV (Ue, Umu)
  and 0.78 GeV (Utau), so that interpolation would lie below the panels and
  nothing is drawn there. Against v1 the upper edge sits 3.5–7× higher for
  Ue and Umu at 0.3–0.9 GeV (the meson-driven `p ↔ n` conversion caps the
  lifetime at ~0.02 s) and 1.7–2.7× lower for Utau (no charged-current
  decay to a charged meson below `m_τ + m_π`, and energetic neutrinos
  neglected, conservatively by the authors' account). Cite Chen and Zhang
  for this layer.
- **The y floor is `1e-10`.** The BBN calculation stops at 1.2 GeV with a
  dashed closure; at `1e-10` the region's upper edge crosses the frame
  floor before that (m = 1.01 / 1.01 / 0.78 GeV for Ue/Umu/Utau), so no
  vertical cliff is drawn. SHiP's minimum in the window, `1.27e-10`, grazes
  the axis on BC6 and BC7; it is the genuine minimum, not a truncation.
- **FASER2's top closure is dropped.** The HNLimits-bundled FASER2
  contours are closed along the top of their source frame (a flat run at
  `1e-3` for Ue/Umu, `1e-1` for Utau, m = 0.1–1.8 GeV). That is a frame
  edge, not a sensitivity boundary; `_drop_top_closure` splits the contour
  at its flat maximum and omits it.

**SHiP and FASER2 are no longer the HNLimits curves (2026-10-25).** The
HNLimits SHiP curves (arXiv:1811.00930) are for 2e20 protons on target, a
third of the exposure of every other SHiP curve in the paper, and its
FASER2 curves (arXiv:1811.12522) use the 2018 strawman geometry without
kaon production. The paper panels now draw:

- SHiP at 6e20: the ESPP 2026 Briefing Book (arXiv:2511.03883v1, Fig. 8.19)
  for Ue and Umu, the PBC ECN3 report (arXiv:2310.17726v1, BC8.pdf) for Utau
  (`curves.tools.vector_ship_hnl`, `data/hnl/vector/SHiP_bb_*.dat`,
  `SHiP_ecn3_Utau.dat`). They reach 1.1-3x smaller couplings but lower masses
  (5.1 GeV for e and mu, 3.5 GeV for tau), consistent with no B_c production.
- FASER2: Feng, Hewitt, Kling, La Rocco, arXiv:2405.07330v2 Fig. 9 (FPF
  geometry, 650 m, 3 x 1 x 10 m; kaon production included), digitised from
  its bitmaps (`curves.tools.raster_faser2_hnl`,
  `data/hnl/raster/FASER2_2405_*.dat`). Up to 10x deeper below 0.6 GeV,
  similar above, mass reach 3.25-3.3 GeV.

The HNLimits files stay for the diagnostic layout. Note that HNLimits 1.2.0
`SHiP/Ue4.dat` copies the muon curve over 0.72-1.85 and 2.80-5.78 GeV
(`Ue4_top.dat` equals `Umu4_top.dat`); the effect is <= 8 %.

CODEX-b is absent from HNLimits, and the PBC BC7 compilation redraws it for
muons only, so its projection is vector-extracted from the EoI's own
per-flavour panels (`curves.tools.vector_codexb_hnl`). The EoI quotes
Dirac HNLs. With equal production and `Γ_M = 2 Γ_D` at fixed `U^2`,
`N_M(U^2) = N_D(2U^2)/2`, so the Majorana `N0` contour is the Dirac `2 N0`
contour at half the coupling. The stored curve follows the CODEX-b
authors' prescription (S. Knapen, 2026-10-05): ×1/√2 on the lower edge,
exact there since `N ∝ U^4`, and ×1/2 on the lifetime-limited upper edge,
the tip staying at the Dirac tip mass; where the upper edge is above the EoI
frame the curve is left open. ×1/2 keeps the lifetime but drops the halved
production, so it is the outer bound of the Majorana island: the exact upper
edge lies lower by a relative `~ln2/(a U^2)` (`a` the exponential slope of
the yield), a few per cent. From 2026-09-28 to 2026-10-06 the curve modelled
that shift per mass column with `N_D = C x^2 exp(-a x)` through the two
Dirac edges, 3–14 % below the bound on the e and μ upper edges and closing
the τ island at 3.07 GeV. v1 applied HNLimits' uniform `dirac_to_majorana_dic["BD"] = 1/√2` to
every vertex, the frame vertex at `U^2 = 1e-2` included (the curve stopped
mid-panel at `7.07e-3`). The unconverted value is the third column. Best
reach as drawn in the EoI (Dirac): `1.88e-8 / 1.99e-8 / 8.38e-7` for
Ue/Umu/Utau, matching the EoI's stated ~1e-8 (e, μ) and ~1e-6 (τ);
Majorana `1.33e-8 / 1.40e-8 / 5.93e-7`. The caption should call the CODEX-b
curve a Majorana recast of the EoI's Dirac projection.

The run conditions ride as the legend title rather than a figure title:
they are a property of the GRENDEL curves alone (SHiP is an SPS beam dump,
FASER2 sits in the forward region, CODEX-b is a transverse LHC detector).
The legend does not repeat the arXiv references, so the caption must cite
SHiP (arXiv:2511.03883; arXiv:2310.17726 for tau), FASER2 (arXiv:2405.07330), CODEX-b
(arXiv:1911.00481) and Chen and Zhang (arXiv:2410.07343, the BBN layer; v1
cited Sabti et al.), state the background-free assumption,
and name the panels (BC6 = `U_eN` top left, BC7 = `U_μN` top right,
BC8 = `U_τN` bottom left).

## BC4 panel

- The closed island begins at the dimuon threshold: `grendel_m_min = 0.211`
  drops the open-topped `ee`-only rows below `2 m_μ`, and the view starts
  at 0.25 GeV to trim the very narrow tip where the upper edge overshoots.
  The CSV itself is unfiltered.
- The frame is the Briefing Book's BC4 frame, `1e-14` to `1e-4`: SHiP's
  lower tip (`7.9e-13`) stays in frame and the decade below it holds the
  legend, wrapped to two rows.
- Projection files are split into separate chains at frame-clip crossings
  (`pdfio.clip_chains`): dropping off-frame vertices in place leaves the
  two survivors adjacent and the polyline plotter bridges them with a
  spurious chord.
- SHiP and FASER2 (14 TeV, 3 ab⁻¹) are the curves of the BC4 panel of the
  ESPP 2026 Physics Briefing Book (arXiv:2511.03883v1, Fig. 8.17 top,
  `Br(h → SS) = 0`), the newest compilation of BC4 projections, read from
  the vector page of the arXiv PDF (there is no source bundle) with both
  axes calibrated from that panel's own ticks (`vector_bc4.briefing_book`).
  The caption cites the Briefing Book for both. The panel's other curves
  (LHCb, MAPP, PREFACE, NA62, DarkQuest, FOREHUNT) are not dedicated LLP
  experiments at the LHC or SHiP, and it has no CODEX-b, ANUBIS or MATHUSLA
  curve; CODEX-b stays the curve of its Expression of Interest, which the
  ECN3 BC4 figure also draws.
- Two earlier SHiP renderings are kept as references and not drawn: the
  ECN3 report's (arXiv:2310.17726, Fig. 24, red dashed; `ship_bc4.dat`),
  whose lower edge lies 1.2–1.36× below the Briefing Book's at 0.3–1.5 GeV,
  equal at 2 and 4 GeV and 1.8× above it at 3 GeV, and the arXiv source of
  2504.06692 (`ship_bc4_2504_reference.dat`, v1's), which agrees with the
  Briefing Book below 1 GeV.
- The gray region against the Briefing Book's (checked on a 160 × 180 grid
  over the paper frame): the same except a region at 0.285–0.44 GeV,
  `sin²θ` 5e-8–5e-7, gray in the Briefing Book only (the PBC 2025 figure
  leaves it open, with no label box over it; the Briefing Book does not say
  which search it adds), and one-cell shifts of the J/ψ and ψ(2S) veto
  stripes and differences above 4 GeV, where GRENDEL has closed. The PBC
  2025 region, shared with BC5, is kept.
- The FASER2 curve of the PBC BC5 figure is tracked as a BC5 reference and
  never drawn in the paper profile (the talk BC4 panel shows it labelled
  "(BC5)"): BC5's quartic coupling also opens `B → K S S`, so a B-driven
  experiment's BC5 reach is 10–60× deeper at 0.5–2 GeV than its BC4 reach,
  converging by ~4 GeV (see `vector_bc4`; SHiP's BC5 curve is
  `data/bc5/projections/ship_pbc_bc5.dat`).

## BC5 panel

- Same axes and conventions as BC4 over the PBC 2025 frame (0.1 GeV to
  m_h/2, `sin²θ` 1e-18 to 1; the source axis is `sin θ`); the excluded
  region is the BC4 file, since the figure's limits constrain the mixing
  angle only.
- The LHC searches for Higgs decays to long-lived particles are added from
  A. Escalante del Valle, arXiv:2509.02564 (JHEP 12 (2025) 129, Fig. 5),
  the one work that translates them to the BC5 plane at BR(h→SS) = 0.01,
  from the data it publishes (github.com/aescalante/LLP_DarkScalarReinterpretation
  at `5788e478`; `curves.tools.lhc_bc5`, `data/bc5/existing/`): CMS
  muon-system showers (arXiv:2402.01898), CMS displaced dimuons
  (arXiv:2112.13769, arXiv:2402.14491) and the ATLAS/CMS LLP summary for
  S → b b. Each is a band over the source's mass points (0.4–60 GeV) joined
  in log-log; the source gives one interval per search and mass and
  nothing between 60 and 62.5 GeV. No compilation covers BC5 the way
  HNLimits covers the HNL, so the PBC 2025 region and these bands are the
  smallest complete set, and the caption cites both (the ESPP 2026 Briefing
  Book, Fig. 8.17, re-plots the same bands). All exclusions are drawn in
  the translucent gray of the other panels, so overlapping limits darken
  and each stays visible; the CMS dimuon band is split at the J/psi-psi(2S)
  and Upsilon windows its search vetoes.
- The four projections are the ANUBIS, CODEX-b, MATHUSLA40 and SHiP
  contours of the same figure, all at BR(h→SS) = 0.01 (`vector_bc5`):
  FASER2, FLArE, MAPP, MAPP2 and the LHCb downstream recast are not drawn
  (dedicated LLP experiments, and the same compact set on every panel).
  Each is the width-2 stroke of its colour; the width-1 stroke is the
  label's leader line and is dropped. The PBC line styles encode the
  maturity of the background estimate and are recorded in the file
  headers, not reproduced: each experiment keeps the style it has on the
  other panels.
- Curves are stored as published, including SHiP's closure along
  m_S = 0.021 GeV (off frame) and `sin θ` = 0.0316, and ANUBIS's
  16-segment polyline; the renderer's x range crops the off-frame part.
- ANUBIS is drawn from its own paper instead (`vector_anubis_bc5`,
  `anubis_conservative_bc5.dat`): arXiv:2510.26932v2 Fig. 14a, ceiling,
  1+2DVs, 95% CL, conservative background, 10-60 GeV. It supersedes the
  collaboration's PBC 2025 curve (`anubis_pbc_bc5.dat`, kept for
  reference) and matches the conservative 1DV projection of the same paper
  in the Higgs-plane figure. The background-free contour of Fig. 14a is
  written to `anubis_bkgfree_bc5.dat` for reference.
- The legend wraps to four columns (`legend_ncol`); the BC4 and BC10
  panels keep their single row, so their published PNGs are unchanged.

**Log axes are calibrated from tick geometry, never from tick-label text
boxes.** In a superscripted label such as `10⁻⁶` the mantissa glyph sits
~1 pt below the tick it annotates, and a bbox-based zero point biases the
whole curve by a constant factor (6–9 % here) while leaving its shape
perfect. `major_y_calibration` derives the zero point from the tick marks
and asserts even spacing so a different tick pattern fails loudly
(`major_x_calibration` does the same for an x axis with two major ticks,
checking the minor ones instead). `vector_bc10`'s 2501.04525 axes moved to
the major ticks on 2026-10-07 (the label-text y was 0.13 pt, +0.40 % in 1/f, off).

## BC10 panel

- Axis: `g_Y = v_h/f` with `f` in the BNT convention (`g_aff = c_f m_f/f`,
  `c_f = 1`) used throughout the ALP model. The reference curves store
  `1/f_BNT`; the renderer scales by `y_scale = 246 = v_h`. PBC
  arXiv:1901.09966 writes `g_Y = 2 v f_l⁻¹` for `f_GKOZ = 2 f_BNT`, so
  `2v_h/f` is correct only when `f` silently means `f_GKOZ`.
- The `η` and `η′` pole windows are unsupported rows (non-finite `peak_N`)
  and are never bridged; nor are finite insensitive intervals. The paper
  bands them from the computed row below to the computed row above.
- Existing bounds: the GKOZ LHCb recast and the 2501.04525 past beam dumps
  above 0.6 GeV; below, the FIPs compilation (arXiv:2305.01715, polygons
  identical in arXiv:2510.05257) with the NA62 regions arXiv:2510.05257 adds.
  Below 0.6 GeV `a -> mu mu` is >99.5 % of the width, so pre-GKOZ and revised
  phenomenology agree there; the compilation also covers the 0.38-0.51 GeV
  gap of the GKOZ recast.
- Projections: SHiP from 2501.04525 (revised phenomenology, 0.6-2 GeV);
  below the η, where 2501.04525 has none, the SHiP curve of the PBC ECN3
  figure (arXiv:2310.17726, `vector_bc10_ecn3`), within 3 % of 2501.04525 at
  0.62-0.82 GeV, drawn as the same curve. The SHiP (ECN3 LoI) curve of
  arXiv:2510.05257 is 1.41x higher in g_Y (no B -> K^(i) tower) and is not
  used. CODEX-b and FASER2 from arXiv:2510.05257; the CODEX-b curve equals
  arXiv:2505.05952's own figure to <1 % (g_Y = v_h c_f/Λ). ANUBIS and
  MATHUSLA have published no BC10 projection; the ESPP 2026 Briefing Book
  (arXiv:2511.03883) has no BC10 figure.
- FASER2 in the paper panel (2026-10-25): `faser2_senscalc_fpf.dat`, computed
  with SensCalc eb2c456 (the corrected widths of the SHiP curve) for its
  "FASER2-FPF" setup (650 m, 3 m x 1 m x 10 m), all LHC production modes,
  N_ev >= 3, 3 ab^-1 at SensCalc's 13 TeV spectra; run headless in
  `campaigns/bc10/faser2_senscalc_20261025`, whose pipeline reproduces
  `ship_senscalc.dat` within 0.5 %. The compiled curve
  (`faser2_fipslhcb25.dat`, arXiv:1811.12522) has the 2018 geometry and
  keeps only e, mu, tau, c, b decays (its width is 5-186x too small above
  1 GeV); it stays for the talk and diagnostic layouts.
- Numerical checks (2026-10-07): the FIPs 2022 and FIPs@LHCb 2025 CHARM/kaon
  polygons, calibrated independently (frame corners vs ticks), agree to
  1.000; the PBC ECN3 figure draws that boundary 2-3 % lower in g_Y at
  0.3-0.8 GeV. Above 0.6 GeV the 2501.04525 past beam dumps exclude less
  than the older CHARM adaptation of the compilations (up to 2.4x higher
  lower edge at 0.65-0.8 GeV).
- The 2501.04525 layers are clipped to that figure's frame, m_a =
  0.6–2.0 GeV: the files keep the path data drawn beneath the PDF clip,
  which was not visually published. SHiP and DarkQuest are split at the
  clip (`x_clip`), the filled past-beam-dump region is clipped with a
  data-coordinate rectangle (v1 filled it from 0.555 GeV).
- **Narrow-island guard.** The central curve has a two-point sensitive
  island at 1.40–1.41 GeV (`peak_N` 3.19 / 3.16, 6 % above threshold,
  flanked by 2.97 and 2.71) that draws as a detached speck between the 1.25
  and 1.47 GeV lobes; the `N>=10` contour already excludes it. It is
  suppressed with `min_island_GeV = 0.05` (BC10 only) and reported on
  stderr. A presentation choice, not a correction: an independent
  six-million-event control at 1.40 GeV reproduces the pocket, and the CSV
  keeps it. The narrowest genuine BC10 island is 0.27 GeV wide; BC4 is a
  single 3.6 GeV segment the guard cannot reach.

The LHCb layer is an existing-exclusion recast (GKOZ arXiv:2310.03524),
not a projection; no ATLAS, CMS or LHCb projection is drawn on any panel
(dedicated LLP experiments only).

## Two features deliberately not smoothed

- **The HNL step at `Ue 0.485` / `Umu 0.380` GeV.** It coincides with an
  `n_events` drop of exactly 300 000 and looks like a splice seam, but the
  `K → ℓN` kinematic limits (0.4932 and 0.3880 GeV) fall exactly in those
  grid gaps: it is the kaon production channel closing, worth +145 % (Ue)
  and +1307 % (Umu) in the lower edge. Utau has no such step because
  `m_K < m_τ`.
- **The BC4 island tip at 3.80 GeV** (`peak_N` 3.44). A meson-only control
  rejects that mass, but the published central adds the `Λ_b` pool, which
  legitimately extends closure there. Deleting the point does not trim the
  tip, it truncates the island into a blunt open edge, because that grid
  point *is* the closure.
