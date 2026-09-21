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
1. SHiP -- published in arXiv:2310.17726 (Fig. 9-PhysicsPotential/BC4.pdf,
   red curve), the citable SHiP BC4 (theta^2 vs m_S) reach. Traced from
   the identical, cleaner two-colour rendering in the SHiP source figure
   arXiv:2504.06692 ``BC4.pdf`` (present only in that paper's arXiv source),
   verified to overlay the published contour across the whole range. Both
   are BC4 = Br(h->SS) = 0; do not confuse with the published 2504.06692
   Fig. 1 scalar panel, which is BC5.
2. CODEX-b -- arXiv:1911.00481 ``B_decay_benchmark.pdf`` (dark-blue
   contours), mixing-only B -> K S production, genuinely BC4. The 300 fb^-1
   (outer) contour is plotted; the 10 fb^-1 one is written for reference.
3. FASER2 and 4. SHiP from the PBC BC5 figure -- NOT BC4-valid, kept as BC5
   references only. The argument "FASER2/SHiP cannot produce Higgs bosons,
   so BC5's BR(h->SS) is irrelevant to them" is wrong: the quartic that
   opens h -> SS also opens B -> K S S, which every B-driven experiment
   sees. Measured on SHiP, the BC5 lower branch is 11x/61x/10x/1.8x/1.00
   deeper than BC4 at m_S = 0.5/1.0/2.0/3.0/4.0 GeV, converging only above
   ~4 GeV where the quartic channel closes; the upper branches agree to
   0.65-1.3x since that edge is lifetime-limited. Neither file may be
   plotted as a BC4 projection.

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
from .sources import BC5_PDF, CODEXB_BDECAY_PDF, SHIP_BC4_PDF, require

PROJ_DIR = BC4_DATA / "projections"
VEC_DIR = BC4_DATA / "vector"
DATE = "2026-07-11"
EXCLUDED_DATE = "2026-07-10"
BC5_SOURCE = "Figures/5-PhysicsReach/BSM_benchmarks/pbc_bc5_gray_final.pdf"
SAMPLES = pdfio.linspace_samples(7)

# PBC BC5 figure: major-tick calibration (PDF points) read from the frame.
# x: 10^0 GeV at 252.41, one decade = 141.205; y: 10^-1 sin(theta) at 60.90,
# one decade = 37.3623 (frame bottom 359.798 = 10^-9).
X_DEC0, X_PER_DEC = 252.4115, 141.2047
Y_DEC_M1, Y_PER_DEC = 60.8999, 37.3623
GRAY = (0.5018997192382812, 0.5019760727882385, 0.5018844604492188)
FASER2_GREEN = (0.0, 0.392, 0.0)
SHIP_ORCHID = (0.855, 0.439, 0.839)
DARKBLUE = (0.0, 0.0, 0.545)


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


def text_spans(page):
    return [(s["text"].strip(), s["bbox"]) for b in page.get_text("dict")["blocks"]
            for line in b.get("lines", []) for s in line["spans"]]


def write(path, chains, to_data, header, clip_x=None):
    pdfio.write_chains(path, pdfio.clip_chains(chains, to_data, clip_x), header,
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

    write(PROJ_DIR / "ship_bc4.dat", pdfio.merge_chains(pdfio.collect_strokes(page, (0.0, 0.0, 1.0), t_samples=SAMPLES)),
          to_data,
          ["SHiP BC4 dark-scalar sensitivity",
           "PUBLISHED SOURCE (cite this): arXiv:2310.17726 (PBC, Post-LS3 Options in",
           "ECN3), Fig 9-PhysicsPotential/BC4.pdf, red SHiP curve.",
           "Traced from the identical, cleaner two-colour rendering in",
           "arXiv:2504.06692 figs/BC4.pdf (SHiP source figure, commented out of the",
           "published PDF); verified to overlay the 2310.17726 SHiP contour.",
           f"Date: {DATE}",
           "Source y axis is theta^2 (= sin^2 theta at these mixings).",
           "Clipped to the figure frame m_S in [0.1, 5.5] GeV."],
          clip_x=(0.1, 5.5))


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
              [f"CODEX-b {lumi} dark-scalar sensitivity, mixing-only (B -> K S) production",
               "Extracted from arXiv:1911.00481 (CODEX-b EOI), figures/B_decay_benchmark.pdf",
               f"Date: {DATE}", "Source y axis is sin^2 theta."])


def pbc_bc5_references(page):
    write(PROJ_DIR / "faser2_pbc_bc5.dat",
          pdfio.merge_chains(pdfio.collect_strokes(page, FASER2_GREEN, t_samples=SAMPLES), tol=1.5), pbc_to_data,
          ["FASER2 dark-scalar sensitivity (dark-green solid curve of the PBC BC5 figure)",
           "Extracted from arXiv:2505.00947v2 / pbc_bc5_gray_final.pdf",
           f"Date: {DATE}",
           "THIS IS A BC5 CURVE AND IS NOT VALID AS A BC4 PROJECTION.",
           "It was originally extracted on the argument that FASER2 cannot produce",
           "Higgs bosons, so the BC5 BR(h->SS)=0.01 component is irrelevant to it and",
           "the curve carries over to BC4. That argument is wrong: the same quartic",
           "coupling also opens B -> K S S, a meson decay a B-driven experiment does",
           "see. Measured on SHiP (BC4 vs BC5, same figure family): the BC5 lower",
           "branch is 11x/61x/10x deeper at m_S = 0.5/1.0/2.0 GeV, converging only",
           "above ~4 GeV. Keep as a BC5 reference only.",
           "Source y axis is sin(theta); values below are squared to sin^2(theta)."])
    write(PROJ_DIR / "ship_pbc_bc5_crosscheck.dat",
          pdfio.merge_chains(pdfio.collect_strokes(page, SHIP_ORCHID, t_samples=SAMPLES), tol=1.5), pbc_to_data,
          ["SHiP curve of the PBC BC5 figure -- BC5 reference (plot ship_bc4.dat for BC4)",
           "Extracted from arXiv:2505.00947v2 / pbc_bc5_gray_final.pdf",
           f"Date: {DATE}",
           "This was extracted as a cross-check on ship_bc4.dat, expecting the two to",
           "coincide because SHiP cannot produce Higgs bosons. THE CROSS-CHECK FAILS,",
           "and correctly so: BC5's quartic coupling also opens B -> K S S, which SHiP",
           "does see. Lower-branch ratio BC4/BC5 = 11x/61x/10x/1.8x/1.00 at",
           "m_S = 0.5/1.0/2.0/3.0/4.0 GeV; upper branches agree to 0.65-1.3x since",
           "that edge is lifetime-limited. A large discrepancy here is EXPECTED and is",
           "not evidence of a problem in ship_bc4.dat.",
           "Source y axis is sin(theta); values below are squared to sin^2(theta)."])


def main() -> int:
    ship()
    codexb()
    bc5 = fitz.open(str(require(BC5_PDF)))[0]
    excluded(bc5)
    pbc_bc5_references(bc5)
    return 0


if __name__ == "__main__":
    sys.exit(main())
