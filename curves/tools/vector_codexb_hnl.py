"""CODEX-b HNL projections from the CODEX-b EoI, all three flavours.

The EoI (arXiv:1911.00481v2) ships no machine-readable contour, so the
curve is pulled from the PDF path operators of its three "Combined"
single-flavour panels, ``figures/HNL_{Ele,Mu,Tau}.pdf`` (labels
``fig:{Ne,mu,tau}_combined``). These are the collaboration's own per-flavour
figures rather than a compilation redraw: all three flavours are present,
CODEX-b is a single stroke colour verified against its legend swatch, and
each panel has an uncluttered pair of log axes. The 2025 CODEX-b ESPP
update (arXiv:2505.05952) contains no HNL contours of its own and defers
to this EoI.

Calibration is from the tick marks, never from tick-label text boxes (a
superscripted label's mantissa sits ~1 pt below the tick it annotates and
biases the whole curve by a constant factor). The major ticks are
recovered geometrically and cross-checked against the frame.

The EoI caption quotes DIRAC HNLs; the GRENDEL curve and the HNLimits
bounds are Majorana. Column 2 applies the factor HNLimits uses for a
production x decay ("beam dump") limit, ``dirac_to_majorana_dic["BD"] =
1/sqrt(2)``; column 3 keeps the unconverted value.

Output: ``curves/data/hnl/original_sources/CODEX-b/CODEX-b_2019_{Ue,Umu,Utau}.dat``
(``mass_GeV  U2_majorana  U2_dirac_raw``).

    python -m curves.tools.vector_codexb_hnl
"""
from __future__ import annotations

import sys
from math import sqrt

import fitz
import numpy as np

from ..plot.data import HNL_DATA
from . import pdfio
from .sources import CODEXB_HNL_PDFS, require

DATE = "2026-08-28"
ARXIV = "arXiv:1911.00481v2"
OUT_DIR = HNL_DATA / "original_sources" / "CODEX-b"
SAMPLES = pdfio.linspace_samples(8)

# Confirmed against the legend swatch immediately left of the "CODEX-b" text.
CODEXB_RGB = (1.0, 0.55, 0.0)
DIRAC_TO_MAJORANA = 1.0 / sqrt(2.0)
PANELS = {"Ue": (r"|U_{eN}|^2", "fig:Ne_combined"), "Umu": (r"|U_{\mu N}|^2", "fig:mu_combined"),
          "Utau": (r"|U_{\tau N}|^2", "fig:tau_combined")}
X_LABELS = (0.2, 1.0, 10.0)          # major x ticks, log axis
Y_TOP_DECADE, Y_N_DECADES = -2.0, 8  # major y ticks run 10^-2 (top) .. 10^-10


def _segments(page):
    """Every black straight-line segment on the page (frame + tick marks)."""
    out = []
    for d in page.get_drawings():
        if d["type"] != "s" or not d["color"] or max(d["color"]) > 0.01:
            continue
        for it in d["items"]:
            if it[0] == "l":
                a, b = it[1], it[2]
                out.append((a.x, a.y, b.x, b.y))
    return out


def calibrate(page):
    """(to_data, frame, meta) using tick geometry only."""
    segs = _segments(page)
    xs = [s[0] for s in segs] + [s[2] for s in segs]
    ys = [s[1] for s in segs] + [s[3] for s in segs]
    fl, fr, ft, fb = min(xs), max(xs), min(ys), max(ys)
    # Ticks are short segments springing from a frame edge; the major ones
    # are the longer of the two lengths present.
    y_ticks = [((s[1] + s[3]) / 2, abs(s[2] - s[0])) for s in segs
               if abs(s[1] - s[3]) < 0.02 and 0.5 < abs(s[2] - s[0]) < 12 and abs(min(s[0], s[2]) - fl) < 1.0]
    x_ticks = [((s[0] + s[2]) / 2, abs(s[3] - s[1])) for s in segs
               if abs(s[0] - s[2]) < 0.02 and 0.5 < abs(s[3] - s[1]) < 12 and abs(max(s[1], s[3]) - fb) < 1.0]
    y_major = sorted({v for v, w in y_ticks if w == max(w for _, w in y_ticks)})
    x_major = sorted({v for v, w in x_ticks if w == max(w for _, w in x_ticks)})
    if len(y_major) != Y_N_DECADES + 1:
        raise RuntimeError(f"expected {Y_N_DECADES + 1} major y ticks, got {len(y_major)}")
    if len(x_major) != len(X_LABELS):
        raise RuntimeError(f"expected {len(X_LABELS)} major x ticks, got {len(x_major)}")
    y_per_dec = float(np.mean(np.diff(y_major)))
    y_top = y_major[0]
    spread = np.ptp(np.diff(y_major))
    if spread > 0.1:
        raise RuntimeError(f"y ticks not evenly spaced (spread {spread:.3f} pt)")
    slope, intercept = np.polyfit(np.log10(X_LABELS), x_major, 1)
    resid = np.max(np.abs(np.polyval([slope, intercept], np.log10(X_LABELS)) - x_major))
    if resid > 0.25:
        raise RuntimeError(f"x tick log fit residual {resid:.3f} pt -- axis may not be log")

    def to_data(x, y):
        return 10.0 ** ((x - intercept) / slope), 10.0 ** (Y_TOP_DECADE - (y - y_top) / y_per_dec)

    return to_data, (fl, fr, ft, fb), {
        "x_per_decade": slope, "y_per_decade": y_per_dec,
        "frame_m": (10.0 ** ((fl - intercept) / slope), 10.0 ** ((fr - intercept) / slope))}


def extract(flavor):
    page = fitz.open(str(require(CODEXB_HNL_PDFS[flavor])))[0]
    to_data, _, meta = calibrate(page)
    chains = []
    for d in page.get_drawings():
        if d["type"] != "s" or not d["color"] or not pdfio.color_matches(d["color"], CODEXB_RGB):
            continue
        pts = pdfio.flatten(d["items"], SAMPLES)
        if len(pts) > 5:            # the 2-point legend swatch is not the curve
            chains.append(pts)
    if len(chains) != 1:
        raise RuntimeError(f"{flavor}: expected 1 CODEX-b chain, found {len(chains)}")
    rows = [to_data(x, y) for x, y in chains[0]]
    meta["n_vertices"] = len(rows)
    meta["min_u2_dirac"] = min(u for _, u in rows)
    meta["m_at_min"] = min(rows, key=lambda r: r[1])[0]
    return rows, meta


def write(flavor, rows, meta):
    ylab, figlbl = PANELS[flavor]
    fname = CODEXB_HNL_PDFS[flavor].rsplit("/", 1)[1]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / f"CODEX-b_2019_{flavor}.dat"
    lo, hi = meta["frame_m"]
    with path.open("w") as fh:
        fh.write(f"# CODEX-b projected HNL sensitivity, {ylab} single-flavour dominance\n")
        fh.write(f"# Source: {ARXIV} (CODEX-b Expression of Interest), {figlbl}\n")
        fh.write(f"# Source figure PDF: sources/arXiv-1911.00481/{fname}\n")
        fh.write("# Method: vector extraction from Mathematica PDF path operators;\n")
        fh.write("#   log-axis calibration from TICK GEOMETRY (not tick-label text bboxes).\n")
        fh.write(f"#   x {meta['x_per_decade']:.4f} pt/decade, y {meta['y_per_decade']:.4f} pt/decade;\n")
        fh.write(f"#   frame spans m = {lo:.4g} .. {hi:.4g} GeV.\n")
        fh.write("# Curve colour: RGB(1.00,0.55,0.00), confirmed against the legend swatch.\n")
        fh.write("# Nature: the EoI quotes DIRAC HNLs. Column 2 applies the Majorana\n")
        fh.write("#   conversion HNLimits uses for production x decay limits,\n")
        fh.write(f"#   dirac_to_majorana_dic['BD'] = 1/sqrt(2) = {DIRAC_TO_MAJORANA:.6f}.\n")
        fh.write("#   Column 3 is the unconverted value as drawn in the source figure.\n")
        fh.write("# Contour is open: it exits the frame at the left edge (low mass) and\n")
        fh.write("#   through the top edge (short lifetime); the region is closed by the\n")
        fh.write("#   plot border in the original, so no closure chord is stored here.\n")
        fh.write(f"# Best reach (Dirac, as drawn): U^2 = {meta['min_u2_dirac']:.4e} at m_N = {meta['m_at_min']:.4g} GeV\n")
        fh.write(f"# Extracted: {DATE}\n")
        fh.write("# mass_GeV  U2_majorana  U2_dirac_raw\n")
        for m, u2 in rows:
            fh.write(f"{m:.8e}  {u2 * DIRAC_TO_MAJORANA:.8e}  {u2:.8e}\n")
    print(f"wrote {path}  ({len(rows)} vertices)")


def main() -> int:
    for flavor in PANELS:
        rows, meta = extract(flavor)
        write(flavor, rows, meta)
        print(f"   {flavor}: {meta['n_vertices']} pts, best Dirac U^2={meta['min_u2_dirac']:.3e} at "
              f"{meta['m_at_min']:.3g} GeV -> Majorana {meta['min_u2_dirac'] * DIRAC_TO_MAJORANA:.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
