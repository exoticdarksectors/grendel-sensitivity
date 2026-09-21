# Figure provenance and conventions

The three sensitivity figures of arXiv:2609.00152 are rendered by
`python -m curves.plot paper` from a results tree holding the central
`N_sig >= 3` scan and the same scan re-solved at `N_sig >= 10`:

| Paper figure | File | Layers |
|---|---|---|
| Fig. 4 | `bc4_grendel_paper.png` | current exclusion (PBC 2025 BC5 figure), SHiP, CODEX-b 300 fb⁻¹, GRENDEL |
| Fig. 5 | `hnlimits_grendel_paper.png` | HNLimits experimental and cosmological bounds, SHiP, FASER2, CODEX-b, GRENDEL; one panel per flavour |
| Fig. 6 | `bc10_grendel_paper.png` | LHCb `B → K μμ` recast (GKOZ), past beam dumps, SHiP, GRENDEL |

The published PNGs have SHA-256

```
bc4_grendel_paper.png       7da7e47e5c741c0580f398f00d4dcfbc9cb058928c29e8bc0ea07424e7e9ec21
hnlimits_grendel_paper.png  dcad42ae1fc6aa206a289b8f27bae6368b8fb86af94aa593fdc6d41f6d71e71a
bc10_grendel_paper.png      d73ee0f8423e0a6937600c203d6e6b66c0ac899520beced66bc750fafcce68cd
```

and the renderer reproduces them byte for byte from the published curves
on Python 3.11 with matplotlib 3.10.9, numpy 2.4.6 and HNLimits 1.2.0.
The PNG is the deterministic artefact; matplotlib stamps a creation date
into the PDF written alongside it, so PDF hashes change on every run.

The paper profiles select layers and the final-size layout only. They
call the same loaders, metadata transformations, clipping rules and
central-curve renderer as the talk and diagnostic layouts, so the two
cannot diverge.

## Two signal thresholds

Every panel draws the `N_sig >= 3` island as the filled red band and the
`N_sig >= 10` island as a dashed contour nested inside it. The panels and
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
(`N ~ (1/f)⁴`) predicts `(10/3)^(1/4) = 1.351`; the published curves give
1.826 and 1.353 at their minima.

Island tips are closed by interpolating the yield-threshold crossing in
log yield and log coupling between adjacent grid points
(`insert_threshold_tips`, `append_high_mass_tip`); a curve solved at a
different threshold closes at that threshold.

## HNL panels (BC6–BC8)

The experimental and cosmological bounds are HNLimits' curves transformed
with the unit, confidence-level and Dirac/Majorana metadata of its
workbook (`curves/data/hnl/local_HNL_database.xlsx`, SHA-256
`aa9346350fd2354af29b71c3a4c6d2da1c5454abd6cef3b780dd1dd858f56783`); the
raw `.dat` files alone lose those transformations. The transformed
regions are written out by `curves.tools.hnlimits_process` and checked in
CI against the tracked copies.

Three local decisions depart from what HNLimits would draw:

- **Cosmology is per flavour.** HNLimits' `cosmo` row points every
  scenario sheet at the same pair of files, `Sabti_BBN/Ue4_{top,bottom}_Cosmo.dat`,
  whose `_bottom` holds two points at one value (the terminal value of the
  top curve, a flat closing line rather than a boundary). Sabti et al.
  (arXiv:2006.07387) publish per-flavour bounds and HNLimits ships them
  alongside as `{Ue4,Umu4,Utau4}_{top,bottom}.dat`; `COSMO_PER_FLAVOUR_FILES`
  substitutes them. `Ue4_top_Cosmo.dat` equals `Ue4_top.dat` except below
  1 MeV, the CMB component exists only below 100 MeV, and the paper quotes
  its bounds for two degenerate Majorana HNLs, so HNLimits' Dirac tag and
  its factor of 1/2 remain correct. Cite Sabti et al. for this layer.
- **The y floor is `1e-10`.** The Sabti calculation stops near the QCD
  scale; with a lower floor the region ended mid-panel in a vertical cliff
  that reads as a physical boundary. At `1e-10` its upper edge crosses the
  frame floor before the data runs out (m = 0.70 / 0.71 / 0.91 GeV for
  Ue/Umu/Utau). SHiP's minimum in the window, `1.27e-10`, grazes the axis
  on BC6 and BC7; it is the genuine minimum, not a truncation.
- **FASER2's top closure is dropped.** The HNLimits-bundled FASER2
  contours are closed along the top of their source frame (a flat run at
  `1e-3` for Ue/Umu, `1e-1` for Utau, m = 0.1–1.8 GeV). That is a frame
  edge, not a sensitivity boundary; `_drop_top_closure` splits the contour
  at its flat maximum and omits it.

CODEX-b is absent from HNLimits and from the PBC BC7 compilation (which is
muon-only), so its projection is vector-extracted from the EoI's own
per-flavour panels (`curves.tools.vector_codexb_hnl`). The EoI quotes
Dirac HNLs; the stored curve applies HNLimits' production × decay factor
`dirac_to_majorana_dic["BD"] = 1/√2`, keeping the unconverted value in a
third column. Best reach as drawn: `1.88e-8 / 1.99e-8 / 8.38e-7` for
Ue/Umu/Utau, matching the EoI's stated ~1e-8 (e, μ) and ~1e-6 (τ).

The run conditions ride as the legend title rather than a figure title:
they are a property of the GRENDEL curves alone (SHiP is an SPS beam dump,
FASER2 sits in the forward region, CODEX-b is a transverse LHC detector).
The legend does not repeat the arXiv references, so the caption must cite
SHiP (arXiv:1811.00930), FASER2 (arXiv:1811.12522), CODEX-b
(arXiv:1911.00481) and Sabti et al., state the background-free assumption,
and name the panels (BC6 = `U_eN` top left, BC7 = `U_μN` top right,
BC8 = `U_τN` bottom left).

## BC4 panel

- The closed island begins at the dimuon threshold: `grendel_m_min = 0.211`
  drops the open-topped `ee`-only rows below `2 m_μ`, and the view starts
  at 0.25 GeV to trim the very narrow tip where the upper edge overshoots.
  The CSV itself is unfiltered.
- The y floor `1e-13` keeps SHiP's lower tip (`7.9e-13`) in frame.
- Projection files are split into separate chains at frame-clip crossings
  (`pdfio.clip_chains`): dropping off-frame vertices in place leaves the
  two survivors adjacent and the polyline plotter bridges them with a
  spurious chord.
- SHiP is cited as arXiv:2310.17726 (the PBC ECN3 report, where the curve
  is published); it is traced from the identical, cleaner rendering in the
  arXiv source of 2504.06692, verified to overlay the published contour.
- The FASER2 and SHiP curves of the PBC BC5 figure are tracked as BC5
  references only and never drawn as BC4 projections: BC5's quartic
  coupling also opens `B → K S S`, so a B-driven experiment's BC5 reach is
  10–60× deeper below ~4 GeV than its BC4 reach (see `vector_bc4`).

**Log axes are calibrated from tick geometry, never from tick-label text
boxes.** In a superscripted label such as `10⁻⁶` the mantissa glyph sits
~1 pt below the tick it annotates, and a bbox-based zero point biases the
whole curve by a constant factor (6–9 % here) while leaving its shape
perfect. `major_y_calibration` derives the zero point from the tick marks
and asserts even spacing so a different tick pattern fails loudly. The
remaining bbox-calibrated extractor is `vector_bc10`'s x axis (0.13 pt =
+0.40 % in 1/f), invisible at plot scale; fix it if that curve is ever
regenerated.

## BC10 panel

- Axis: `g_Y = v_h/f` with `f` in the BNT convention (`g_aff = c_f m_f/f`,
  `c_f = 1`) used throughout the ALP model. The reference curves store
  `1/f_BNT`; the renderer scales by `y_scale = 246 = v_h`. PBC
  arXiv:1901.09966 writes `g_Y = 2 v f_l⁻¹` for `f_GKOZ = 2 f_BNT`, so
  `2v_h/f` is correct only when `f` silently means `f_GKOZ`.
- The `η` and `η′` pole windows are unsupported rows (non-finite `peak_N`)
  and are never bridged; nor are finite insensitive intervals.
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
