"""BC4 (dark scalar, mixing only) reference curves.

Existing exclusion
------------------
The 2025 PBC report (arXiv:2505.00947) has no BC4 figure -- that update
covers BC3/BC5/BC7, the PBC BC4 sensitivity figure is in the ECN3 report
arXiv:2310.17726 -- but its BC5 figure (``pbc_bc5_gray_final.pdf``) carries
the dark-scalar existing limits (E949, NA62, MicroBooNE, KOTO, ICARUS, LHCb,
Belle II, plus reinterpretations of PS191, CHARM, LSND). Those constrain
the mixing angle only, so the gray region applies to BC4 unchanged. The
figure's y axis is sin(theta); vertices are squared on output. Calibration
from the figure's own major ticks.

Projections (dedicated LLP experiments only)
--------------------------------------------
1. SHiP and FASER2 -- arXiv:2511.03883 (ESPP 2026 Physics Briefing Book),
   Fig. 8.17 top panel (Br(h->SS) = 0), the newest compilation of BC4
   projections; drawn. The Briefing Book is a PDF without a source bundle,
   so the curves are read from the vector page (``briefing_book``). The
   panel's other projections (LHCb, MAPP, PREFACE, NA62, DarkQuest,
   FOREHUNT) are not dedicated LLP experiments at the LHC or SHiP.
2. SHiP, earlier renderings -- references, never drawn:
   arXiv:2310.17726 (PBC ECN3 report), Fig. 24 =
   ``Figures/9-PhysicsPotential/BC4.pdf``, red dashed curve
   (``ship_bc4.dat``), and the two-colour rendering in the SHiP source
   figure arXiv:2504.06692 ``BC4.pdf`` (present only in that paper's arXiv
   source, not citable; ``ship_bc4_2504_reference.dat``). Against the ECN3
   curve the Briefing Book's lower edge is 1.2-1.36x higher below 1.5 GeV,
   equal at 2 and 4 GeV and 1.8x lower at 3 GeV; the 2504 trace agrees with
   the Briefing Book below 1 GeV. All three are BC4 = Br(h->SS) = 0; do not
   confuse with the published 2504.06692 Fig. 1 scalar panel, which is BC5.
3. CODEX-b -- arXiv:1911.00481 ``B_decay_benchmark.pdf`` (dark-blue
   contours), mixing-only b -> s (B -> X_s S) production, genuinely BC4. The 300 fb^-1
   (outer) contour is plotted; the 10 fb^-1 one is written for reference.
   The ECN3 BC4 figure draws it as its representative large-angle LHC
   experiment; the Briefing Book panel has no CODEX-b, ANUBIS or MATHUSLA
   curve for BC4.
4. FASER2 from the PBC BC5 figure -- NOT BC4-valid, kept as a BC5
   reference only. The argument "FASER2 cannot produce Higgs bosons, so
   BC5's BR(h->SS) is irrelevant to it" is wrong: the quartic that opens
   h -> SS also opens B -> K S S, which every B-driven experiment sees.
   Measured on SHiP (``curves/data/bc5/projections/ship_pbc_bc5.dat`` vs
   ``ship_bc4.dat``), the BC5 lower branch is 11x/61x/10x/1.8x/1.00 deeper
   than BC4 at m_S = 0.5/1.0/2.0/3.0/4.0 GeV, converging only above ~4 GeV
   where the quartic channel closes; the upper branches agree to 0.65-1.3x
   since that edge is lifetime-limited. It is not a BC4 projection (the
   talk BC4 panel shows it labelled "(BC5)"; the paper profile does not
   draw it). Unverified counterpoint (2026-10-25): the ESPP 2026 Briefing
   Book draws the identical FASER2 path in its BC4 and BC5 panels, and the
   curve has no B -> X_s S S feature, so it may be mixing-only after all. The BC5 projections themselves (ANUBIS, CODEX-b, MATHUSLA40,
   SHiP) are ``vector_bc5``.

Log axes are calibrated from TICK GEOMETRY, never from tick-label text
boxes: in a superscripted label such as ``10^-6`` the mantissa glyph is
centred below the tick it annotates, and a bbox-based zero point biases
every value by a constant factor while leaving the shape perfect.

    python -m curves.tools.vector_bc4
"""
from __future__ import annotations

import sys

import fitz
import numpy as np

from ..plot.data import BC4_DATA
from . import pdfio
from .sources import BB_PDF, BB_SCALAR_PAGE, BC5_PDF, CODEXB_BDECAY_PDF, ECN3_BC4_PDF, SHIP_BC4_PDF, require

PROJ_DIR = BC4_DATA / "projections"
VEC_DIR = BC4_DATA / "vector"
DATE = "2026-07-11"
ECN3_DATE = "2026-09-28"
EXCLUDED_DATE = "2026-07-10"
FASER2_DATE = "2026-09-26"     # re-extracted without the label leader line
BB_DATE = "2026-10-09"
BC5_SOURCE = "Figures/5-PhysicsReach/BSM_benchmarks/pbc_bc5_gray_final.pdf"
SAMPLES = pdfio.linspace_samples(7)

# PBC BC5 figure: major-tick calibration (PDF points) read from the frame.
# x: 10^0 GeV at 252.41, one decade = 141.205; y: 10^-1 sin(theta) at 60.90,
# one decade = 37.3623 (frame bottom 359.798 = 10^-9).
X_DEC0, X_PER_DEC = 252.4115, 141.2047
Y_DEC_M1, Y_PER_DEC = 60.8999, 37.3623
GRAY = (0.5018997192382812, 0.5019760727882385, 0.5018844604492188)
FASER2_GREEN = (0.0, 0.392, 0.0)
DARKBLUE = (0.0, 0.0, 0.545)
ECN3_SHIP_RED = (1.0, 0.0, 0.0)
BB_SHIP_BLUE = (0.255, 0.412, 0.882)
BB_FASER2_MAGENTA = (1.0, 0.0, 1.0)


def pbc_to_data(x, y):
    log_m = (x - X_DEC0) / X_PER_DEC
    log_sinth = -1.0 - (y - Y_DEC_M1) / Y_PER_DEC
    return 10 ** log_m, 10 ** (2.0 * log_sinth)


def major_y_calibration(page, decades_between_major):
    """(y of the top major tick, points per decade) from tick geometry.

    Ticks are short marks anchored on the frame's left edge; they are
    stroked line items in some figures and filled rects in others, so both
    are collected. Major ticks are the longer of the two lengths present,
    and their spacing is asserted even so a figure with a different tick
    pattern fails loudly instead of drifting.
    """
    marks = []
    verticals = []
    for d in page.get_drawings():
        black = d.get("color") and max(d["color"]) < 0.01
        if d["type"] == "s" and black:
            for it in d["items"]:
                if it[0] != "l":
                    continue
                a, b = it[1], it[2]
                if abs(a.x - b.x) < 0.05 and abs(a.y - b.y) > 50:
                    verticals.append(a.x)
                elif abs(a.y - b.y) < 0.05 and 0.5 < abs(a.x - b.x) < 25:
                    marks.append(((a.y + b.y) / 2, abs(a.x - b.x), min(a.x, b.x), max(a.x, b.x)))
        if d["type"] in ("f", "fs") and d.get("fill") and max(d["fill"]) < 0.01:
            r = d["rect"]
            if r.height < 0.5 and 0.5 < r.width < 25:
                marks.append(((r.y0 + r.y1) / 2, r.width, r.x0, r.x1))
    if not verticals or not marks:
        raise RuntimeError("could not locate plot frame or y tick marks")
    frame_left = min(verticals)
    at_edge = [m for m in marks if abs(m[2] - frame_left) < 2.0 or abs(m[3] - frame_left) < 2.0]
    if not at_edge:
        raise RuntimeError("no y tick marks anchored on the frame's left edge")
    longest = max(m[1] for m in at_edge)
    major = sorted({round(m[0], 3) for m in at_edge if abs(m[1] - longest) < 0.2})
    if len(major) < 2:
        raise RuntimeError(f"expected >=2 major y ticks, found {len(major)}")
    steps = np.diff(major)
    if float(np.ptp(steps)) > 0.1:
        raise RuntimeError(f"major y ticks unevenly spaced (spread {np.ptp(steps):.3f} pt)")
    return major[0], float(np.mean(steps)) / float(decades_between_major)


def _black_lines(page):
    for d in page.get_drawings():
        if d.get("color") and max(d["color"]) < 0.01:
            yield from ((it[1], it[2]) for it in d["items"] if it[0] == "l")


def frame_edges(page):
    """(left, right, top, bottom) of the plot frame: the long black lines."""
    xs = [a.x for a, b in _black_lines(page) if abs(a.x - b.x) < 0.05 and abs(a.y - b.y) > 50]
    ys = [a.y for a, b in _black_lines(page) if abs(a.y - b.y) < 0.05 and abs(a.x - b.x) > 50]
    return min(xs), max(xs), min(ys), max(ys)


def major_x_calibration(page):
    """(x of the leftmost major tick, points per decade) from tick geometry.

    The x counterpart of ``major_y_calibration`` for an axis with a major
    tick at every decade: short vertical marks springing from the frame's
    bottom edge, the major ones being the longer of the two lengths
    present. Two majors cannot reveal uneven spacing, so every minor tick
    must also land on a log10(k) position (k = 2..9) to within 0.1 pt.
    """
    bottom = frame_edges(page)[3]
    marks = [((a.x + b.x) / 2, abs(a.y - b.y)) for a, b in _black_lines(page)
             if abs(a.x - b.x) < 0.05 and 0.5 < abs(a.y - b.y) < 25 and abs(min(a.y, b.y) - bottom) < 2.0]
    if not marks:
        raise RuntimeError("no x tick marks on the frame's bottom edge")
    longest = max(w for _, w in marks)
    major = sorted({round(x, 3) for x, w in marks if abs(w - longest) < 0.2})
    minor = sorted({round(x, 3) for x, w in marks if abs(w - longest) >= 0.2})
    if len(major) < 2:
        raise RuntimeError(f"expected >=2 major x ticks, found {len(major)}")
    steps = np.diff(major)
    if float(np.ptp(steps)) > 0.1:
        raise RuntimeError(f"major x ticks unevenly spaced (spread {np.ptp(steps):.3f} pt)")
    per_dec = float(np.mean(steps))
    frac = ((np.asarray(minor) - major[0]) / per_dec) % 1.0
    resid = np.min(np.abs(frac[:, None] - np.log10(np.arange(2, 10))[None, :]), axis=1) * per_dec
    if minor and float(resid.max()) > 0.1:
        raise RuntimeError(f"minor x ticks off the log10(k) grid by up to {resid.max():.3f} pt")
    return major[0], per_dec


def text_spans(page):
    return [(s["text"].strip(), s["bbox"]) for b in page.get_text("dict")["blocks"]
            for line in b.get("lines", []) for s in line["spans"]]


def write(path, chains, to_data, header, clip_x=None, clip_y=None):
    pdfio.write_chains(path, pdfio.clip_chains(chains, to_data, clip_x, clip_y), header,
                       columns='Format: chains separated by blank line; each row "mass_GeV  sin2theta"')


def excluded(page):
    main = next(d for d in page.get_drawings() if d.get("type") in ("f", "fs") and d.get("fill")
                and pdfio.near(d["fill"], GRAY) and len(d.get("items", [])) > 100)
    subs = pdfio.subpaths(main["items"])
    print(f"{len(subs)} subpaths, sizes: {[len(s) for s in subs[:10]]} ...")
    VEC_DIR.mkdir(parents=True, exist_ok=True)
    path = VEC_DIR / "currently_excluded.dat"
    total = 0
    with path.open("w") as fh:
        fh.write('# "Currently excluded" gray region from the BC5 fig (mixing-only limits, valid for BC4)\n')
        fh.write(f"# Extracted from arXiv:2505.00947v2 / {BC5_SOURCE}\n")
        fh.write(f"# Date: {EXCLUDED_DATE}\n")
        fh.write("# Source y axis is sin(theta); values below are squared to sin^2(theta).\n")
        fh.write('# Format: subpaths separated by blank line; each row is "mass_GeV  sin2theta"\n')
        fh.write("# Source experiments per caption: E949, NA62, MicroBooNE, KOTO, ICARUS,\n")
        fh.write("# LHCb, BelleII, plus reinterpretations of PS191, CHARM, LSND.\n")
        fh.write("# mass_GeV  sin2theta\n")
        for k, sub in enumerate(subs):
            if k > 0:
                fh.write("\n")
            fh.write(f"# subpath {k}: {len(sub)} vertices\n")
            for x, y in sub:
                m, s2t = pbc_to_data(x, y)
                if m < 1e-2 or m > 1e2 or s2t < 1e-19:
                    continue
                fh.write(f"{m:.5e}  {s2t:.5e}\n")
                total += 1
    print(f"wrote {path} ({total} vertices total)")


def ship():
    """The SHiP curve of the ECN3 report's own figure (2310.17726 Fig. 24), a reference."""
    page = fitz.open(str(require(ECN3_BC4_PDF)))[0]
    # A major tick at every decade on both axes: sin^2 theta 10^-5 (top) ..
    # 10^-12, m_S 0.1 and 1 GeV.
    cy5, ypd = major_y_calibration(page, decades_between_major=1.0)
    cx01, xpd = major_x_calibration(page)

    def to_data(x, y):
        return 10 ** (-1.0 + (x - cx01) / xpd), 10 ** (-5.0 - (y - cy5) / ypd)

    left, right, top, bottom = frame_edges(page)
    (m_lo, u_top), (m_hi, u_bottom) = to_data(left, top), to_data(right, bottom)
    # One dashed stroke holding two subpaths (the upper edge, then the lower
    # edge leftwards from the 5 GeV grid edge): split at the pen jump, rejoin
    # at the shared vertex. The path continues beneath the PDF clip, so it is
    # clipped to the frame on both axes.
    chains = [sub for d in page.get_drawings()
              if d["type"] == "s" and d["color"] and pdfio.color_matches(d["color"], ECN3_SHIP_RED)
              and pdfio.is_dashed(d) for sub in pdfio.subpaths(d["items"], SAMPLES)]
    write(PROJ_DIR / "ship_bc4.dat", pdfio.merge_chains(chains), to_data,
          ["SHiP BC4 dark-scalar sensitivity",
           "Extracted from arXiv:2310.17726 (PBC, Post-LS3 Options in ECN3), Fig. 24 =",
           "Figures/9-PhysicsPotential/BC4.pdf, red dashed SHiP curve. REFERENCE: the paper",
           "draws the Briefing Book curve, ship_bb_bc4.dat.",
           f"Date: {ECN3_DATE}",
           f"Calibration from tick geometry on both axes: x {xpd:.4f} pt/decade (major",
           "ticks 0.1 and 1 GeV, every minor tick on a log10(k) position), y",
           f"{ypd:.4f} pt/decade (major ticks at every decade, 10^-5 at the frame top).",
           "Source y axis is sin^2 theta.",
           f"Clipped to the figure frame m_S in [{m_lo:.3g}, {m_hi:.3g}] GeV, sin^2 theta in",
           f"[{u_bottom:.3g}, {u_top:.3g}]; the vertical edge at m_S = 5.0 GeV is the source's own.",
           "Replaces the trace of the arXiv:2504.06692 source figure, kept as",
           "ship_bc4_2504_reference.dat: that lower edge is 1.39/1.36/1.20x higher at",
           "m_S = 0.5/1/2 GeV; the two agree above ~3 GeV and on the upper edge."],
          clip_x=(m_lo, m_hi), clip_y=(u_bottom, u_top))


def ship_2504_reference():
    """The SHiP source figure of arXiv:2504.06692 -- a reference, never drawn."""
    page = fitz.open(str(require(SHIP_BC4_PDF)))[0]
    spans = text_spans(page)

    def cx(txt):
        for t, bb in spans:
            if t == txt:
                return (bb[0] + bb[2]) / 2
        raise KeyError(txt)

    # x is safe from label bboxes (horizontal centring is exact here), y is
    # not. Major y ticks are 3 decades apart: 10^-6 (top), 10^-9, 10^-12.
    cx01, cx1 = cx("0.1"), cx("1")
    xpd = cx1 - cx01
    cy6, ypd = major_y_calibration(page, decades_between_major=3.0)

    def to_data(x, y):
        return 10 ** ((x - cx1) / xpd), 10 ** (-6.0 - (y - cy6) / ypd)

    write(PROJ_DIR / "ship_bc4_2504_reference.dat",
          pdfio.merge_chains(pdfio.collect_strokes(page, (0.0, 0.0, 1.0), t_samples=SAMPLES)), to_data,
          ["SHiP BC4 dark-scalar sensitivity -- REFERENCE ONLY, not drawn (plot ship_bb_bc4.dat)",
           "Traced from the two-colour rendering in arXiv:2504.06692 figs/BC4.pdf (SHiP",
           "source figure, commented out of the published PDF, so not citable).",
           "Against the published contour (arXiv:2310.17726 Fig. 24, ship_bc4.dat) it",
           "agrees above ~3 GeV and on the upper edge, but its lower edge is",
           "1.39/1.36/1.20x higher at m_S = 0.5/1/2 GeV.",
           f"Date: {DATE}",
           "Source y axis is theta^2 (= sin^2 theta at these mixings).",
           "Clipped to the figure frame m_S in [0.1, 5.5] GeV."],
          clip_x=(0.1, 5.5))


def top_panel_frame(page):
    """(left, right, top, bottom) of the upper of the stacked plot frames on a page.

    Frames are the long black lines; the BC4 panel of the Briefing Book page
    sits above the BC5 one and shares no edge with it.
    """
    hs = sorted({(round(a.y, 2), round(min(a.x, b.x), 2), round(max(a.x, b.x), 2))
                 for a, b in _black_lines(page) if abs(a.y - b.y) < 0.05 and abs(a.x - b.x) > 50})
    top = hs[0]
    bottom = next(h for h in hs[1:] if abs(h[1] - top[1]) < 0.1 and abs(h[2] - top[2]) < 0.1)
    return top[1], top[2], top[0], bottom[0]


def panel_tick_calibration(page, frame):
    """Major-tick positions and points per decade on both axes of one panel.

    The tick marks are short black lines on the frame's left (y) and bottom
    (x) edges, the major ones twice as long as the minor ones; there is a
    major tick at every decade on both axes. Every minor tick must sit on a
    log10(k) position (k = 2..9) to within 0.1 pt, so a misread frame or a
    different tick pattern fails loudly.
    """
    left, right, top, bottom = frame
    xs, ys = [], []
    for a, b in _black_lines(page):
        if abs(a.x - b.x) < 0.05 and 0.5 < abs(a.y - b.y) < 25 and abs(max(a.y, b.y) - bottom) < 1.5 \
                and left - 1 < a.x < right + 1:
            xs.append((a.x, abs(a.y - b.y)))
        if abs(a.y - b.y) < 0.05 and 0.5 < abs(a.x - b.x) < 25 and abs(min(a.x, b.x) - left) < 1.5 \
                and top - 1 < a.y < bottom + 1:
            ys.append((a.y, abs(a.x - b.x)))
    out = []
    for marks in (xs, ys):
        longest = max(w for _, w in marks)
        major = sorted({round(p, 3) for p, w in marks if abs(w - longest) < 0.2})
        minor = sorted({round(p, 3) for p, w in marks if abs(w - longest) >= 0.2})
        steps = np.diff(major)
        if len(major) < 2 or float(np.ptp(steps)) > 0.1:
            raise RuntimeError(f"major ticks missing or unevenly spaced: {major}")
        out.append((major, minor, float(np.mean(steps))))
    (xmaj, xmin, xpd), (ymaj, ymin, ypd) = out
    # y grows downwards: a minor tick lies log10(k) decades above the major tick below it.
    for minor, origin, per_dec, sign in ((xmin, xmaj[0], xpd, 1.0), (ymin, ymaj[-1], ypd, -1.0)):
        frac = ((sign * (np.asarray(minor) - origin)) / per_dec) % 1.0
        resid = np.min(np.abs(frac[:, None] - np.log10(np.arange(2, 10))[None, :]), axis=1) * per_dec
        if float(resid.max()) > 0.1:
            raise RuntimeError(f"minor ticks off the log10(k) grid by up to {resid.max():.3f} pt")
    return xmaj, xpd, ymaj, ypd


def briefing_book():
    """SHiP and FASER2 from the BC4 panel of the ESPP 2026 Physics Briefing Book.

    arXiv:2511.03883v1, Fig. 8.17 (top), "Br(h -> SS) = 0": the newest
    compilation of BC4 projections. The page is vector graphics; both axes
    carry a major tick at every decade (m_S 0.1 and 1 GeV, sin^2 theta 1e-4
    at the frame top down to 1e-14 at its bottom). The curves run beneath
    the frame's clip, so they are clipped to it here.
    """
    page = fitz.open(str(require(BB_PDF)))[BB_SCALAR_PAGE]
    frame = top_panel_frame(page)
    left, right, top, bottom = frame
    xmaj, xpd, ymaj, ypd = panel_tick_calibration(page, frame)
    if len(xmaj) != 2 or len(ymaj) != 11 or abs(xmaj[0] - left) > 0.1 or abs(ymaj[0] - top) > 0.1:
        raise RuntimeError(f"unexpected tick pattern: x {xmaj}, y {ymaj}")

    def to_data(x, y):
        return 10 ** (-1.0 + (x - xmaj[0]) / xpd), 10 ** (-4.0 - (y - ymaj[0]) / ypd)

    (m_lo, u_top), (m_hi, u_bottom) = to_data(left, top), to_data(right, bottom)

    def curve(color):
        # The legend swatches sit right of the frame; the BC5 panel below it.
        chains = [sub for d in page.get_drawings()
                  if d["type"] == "s" and d["color"] and pdfio.color_matches(d["color"], color)
                  and d["rect"].y0 < bottom and d["rect"].x0 < right and d["rect"].width > 20
                  for sub in pdfio.subpaths(d["items"], SAMPLES)]
        return pdfio.merge_chains(chains)

    common = [f"Date: {BB_DATE}",
              f"Calibration from tick geometry: x {xpd:.4f} pt/decade (major ticks 0.1 and 1 GeV),",
              f"y {ypd:.4f} pt/decade (major ticks at every decade, 1e-4 at the frame top);",
              "every minor tick on a log10(k) position. Source y axis is sin^2 theta.",
              f"Clipped to the figure frame m_S in [{m_lo:.3g}, {m_hi:.3g}] GeV, sin^2 theta in",
              f"[{u_bottom:.3g}, {u_top:.3g}]."]
    write(PROJ_DIR / "ship_bb_bc4.dat", curve(BB_SHIP_BLUE), to_data,
          ["SHiP BC4 dark-scalar sensitivity (15 y, 6e20 protons on target)",
           "Extracted from arXiv:2511.03883v1 (ESPP 2026 Physics Briefing Book), Fig. 8.17",
           "top panel (Br(h -> SS) = 0), royal-blue solid SHiP curve (cite this).",
           *common,
           "Against the ECN3 rendering (ship_bc4.dat) the lower edge is 1.2-1.36x higher",
           "at 0.3-1.5 GeV, equal at 2 and 4 GeV and 1.8x lower at 3 GeV; the upper edge",
           "is 1.2-1.4x higher below 1.5 GeV and equal above."],
          clip_x=(m_lo, m_hi), clip_y=(u_bottom, u_top))
    write(PROJ_DIR / "faser2_bc4.dat", curve(BB_FASER2_MAGENTA), to_data,
          ["FASER2 BC4 dark-scalar sensitivity (14 TeV, 3 ab^-1)",
           "Extracted from arXiv:2511.03883v1 (ESPP 2026 Physics Briefing Book), Fig. 8.17",
           "top panel (Br(h -> SS) = 0), magenta solid FASER2 curve (cite this).",
           *common],
          clip_x=(m_lo, m_hi), clip_y=(u_bottom, u_top))


def codexb():
    page = fitz.open(str(require(CODEXB_BDECAY_PDF)))[0]
    # x labels are split into characters; centres of the '0.5' and '1.0'
    # digit groups. y from tick geometry: the '10' spans are doubly unsafe
    # here since the exponent of 10^-10 is itself a "10" span. Major y
    # ticks are 2 decades apart: 10^-6 (top) .. 10^-14.
    cx05 = (120.1 + 132.8) / 2
    cx10 = (174.6 + 187.4) / 2
    xpd = (cx10 - cx05) / 0.3010299957
    cy_m6, ypd = major_y_calibration(page, decades_between_major=2.0)

    def to_data(x, y):
        return 10 ** ((x - cx10) / xpd), 10 ** (-6.0 - (y - cy_m6) / ypd)

    chains = pdfio.merge_chains(pdfio.collect_strokes(page, DARKBLUE, t_samples=SAMPLES), tol=1.0)
    chains.sort(key=lambda c: -len(c))
    # 300 fb^-1 reaches lower sin2theta than 10 fb^-1.
    chains_sorted = sorted(chains, key=lambda ch: min(to_data(x, y)[1] for x, y in ch))
    c300, c10 = chains_sorted[0], chains_sorted[-1]
    for fname, chain, lumi in (("codexb_300fb.dat", c300, "300 fb^-1"), ("codexb_10fb.dat", c10, "10 fb^-1")):
        write(PROJ_DIR / fname, [chain], to_data,
              [f"CODEX-b {lumi} dark-scalar sensitivity, mixing-only (b -> s, B -> X_s S) production",
               "Extracted from arXiv:1911.00481 (CODEX-b EOI), figures/B_decay_benchmark.pdf",
               f"Date: {DATE}", "Source y axis is sin^2 theta."])


def faser2_bc5_reference(page):
    # width=2.0 selects the curve; the width-1 stroke of the same colour is
    # the label's leader line.
    write(PROJ_DIR / "faser2_pbc_bc5.dat",
          pdfio.merge_chains(pdfio.collect_strokes(page, FASER2_GREEN, t_samples=SAMPLES, width=2.0), tol=1.5),
          pbc_to_data,
          ["FASER2 dark-scalar sensitivity (dark-green solid curve of the PBC BC5 figure)",
           "Extracted from arXiv:2505.00947v2 / pbc_bc5_gray_final.pdf",
           f"Date: {FASER2_DATE}",
           "THIS IS A BC5 CURVE AND IS NOT VALID AS A BC4 PROJECTION.",
           "It was originally extracted on the argument that FASER2 cannot produce",
           "Higgs bosons, so the BC5 BR(h->SS)=0.01 component is irrelevant to it and",
           "the curve carries over to BC4. That argument is wrong: the same quartic",
           "coupling also opens B -> K S S, a meson decay a B-driven experiment does",
           "see. Measured on SHiP (BC4 vs BC5, same figure family): the BC5 lower",
           "branch is 11x/61x/10x deeper at m_S = 0.5/1.0/2.0 GeV, converging only",
           "above ~4 GeV. Keep as a BC5 reference only.",
           "Source y axis is sin(theta); values below are squared to sin^2(theta)."])


def main() -> int:
    briefing_book()
    ship()
    ship_2504_reference()
    codexb()
    bc5 = fitz.open(str(require(BC5_PDF)))[0]
    excluded(bc5)
    faser2_bc5_reference(bc5)
    return 0


if __name__ == "__main__":
    sys.exit(main())
