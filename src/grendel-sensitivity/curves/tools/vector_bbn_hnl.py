"""The BBN bound on a single Majorana HNL, all three flavours.

The cosmological layer of the HNL panels. The paper's benchmark is a
single Majorana HNL mixing with one flavour, and the BBN calculations in
common use do not describe it: Sabti et al. (arXiv:2006.07387, the
HNLimits ``cosmo`` row) simulate a Dirac HNL, i.e. two degenerate Majorana
HNLs, and state that "the bounds do not apply to models with single
Majorana HNLs" (their Sec. 5.2); Boyarsky et al. (arXiv:2008.00749) also
take a degenerate pair forming one quasi-Dirac fermion. Chen and Zhang,
"BBN constraint on heavy neutrino production and decay", arXiv:2410.07343
(Phys. Rev. D 111 (2025) 123024), work in exactly the benchmark model:
one Majorana heavy neutrino mixing with a single flavour, whose early-
universe abundance (thermal freeze-out, or freeze-in from zero) and decays
are both fixed by the mass and mixing, excluded at 95% CL. The curve is
used as published: no CL rescaling and no Dirac/Majorana factor.

No machine-readable curve is published, so the blue "BBN Constraint"
region is read from the figure PDFs of the arXiv source: Fig. 1
(``Ue4sq_Constraints.pdf``) and the two panels of Fig. 8
(``Umu4sq_constraints.pdf``, ``Utau4sq_constraints.pdf``). The region is a
30%-opacity navy fill bounded by two solid navy strokes, an upper and a
lower boundary, which are the stored curves; every stroke vertex is checked
to be a vertex of the fill. The fill's dashed right edge at ~1.2 GeV is
where the calculation stops (the paper leaves 1-10 GeV as a gap), not a
boundary. The >= 10 GeV region and, above 0.3 GeV, the CMB constraint lie
at |U|^2 < 1e-15, below every panel, and are not stored.

Both log axes are calibrated from tick geometry: a major y tick at every
decade (10^0 .. 10^-28, labelled every other decade) and major x ticks at
10^-3 .. 10 and 30 GeV. The minor ticks of these figures sit at fifths of
a decade, not at log10(k), so they are not used.

Output: ``curves/data/hnl/cosmology/BBN_ChenZhang2025_{Ue,Umu,Utau}.dat``
(``# bottom`` / ``# top`` sections, ``mass_GeV  U2``).

    python -m curves.tools.vector_bbn_hnl
"""
from __future__ import annotations

import sys

import fitz
import numpy as np

from ..plot.data import HNL_DATA
from . import pdfio
from .sources import BBN_HNL_PDFS, require

DATE = "2026-09-28"
ARXIV = "arXiv:2410.07343v2"
OUT_DIR = HNL_DATA / "cosmology"
NAVY = (0.0, 0.0, 0.502)
PANELS = {"Ue": (r"|U_{eN}|^2", "Fig. 1"), "Umu": (r"|U_{\mu N}|^2", "Fig. 8 (top)"),
          "Utau": (r"|U_{\tau N}|^2", "Fig. 8 (bottom)")}
X_LABELS = (1.0e-3, 1.0e-2, 1.0e-1, 1.0, 10.0, 30.0)   # major x ticks, log axis
Y_TOP_DECADE, Y_N_DECADES = 0.0, 28                    # a major y tick per decade, 10^0 (top) .. 10^-28
LOW_MASS_EDGE_GeV = 5.0                                 # the low-mass region ends at ~1.2 GeV


def _black_lines(page):
    for d in page.get_drawings():
        if d.get("color") and max(d["color"]) < 0.01:
            yield from ((it[1], it[2]) for it in d["items"] if it[0] == "l")


def calibrate(page):
    """(to_data, meta) from the tick marks springing from the frame's left and bottom edges."""
    lines = list(_black_lines(page))
    left = min(a.x for a, b in lines if abs(a.x - b.x) < 0.05 and abs(a.y - b.y) > 50)
    bottom = max(a.y for a, b in lines if abs(a.y - b.y) < 0.05 and abs(a.x - b.x) > 50)
    y_ticks = [((a.y + b.y) / 2, abs(a.x - b.x)) for a, b in lines
               if abs(a.y - b.y) < 0.05 and 0.5 < abs(a.x - b.x) < 12 and abs(min(a.x, b.x) - left) < 0.5]
    x_ticks = [((a.x + b.x) / 2, abs(a.y - b.y)) for a, b in lines
               if abs(a.x - b.x) < 0.05 and 0.5 < abs(a.y - b.y) < 12 and abs(max(a.y, b.y) - bottom) < 0.5]
    y_major = sorted({round(v, 3) for v, w in y_ticks if w == max(w for _, w in y_ticks)})
    x_major = sorted({round(v, 3) for v, w in x_ticks if w == max(w for _, w in x_ticks)})
    if len(y_major) != Y_N_DECADES + 1:
        raise RuntimeError(f"expected {Y_N_DECADES + 1} major y ticks, got {len(y_major)}")
    if len(x_major) != len(X_LABELS):
        raise RuntimeError(f"expected {len(X_LABELS)} major x ticks, got {len(x_major)}")
    spread = float(np.ptp(np.diff(y_major)))
    if spread > 0.1:
        raise RuntimeError(f"y ticks not evenly spaced (spread {spread:.3f} pt)")
    y_per_dec, y_top = float(np.mean(np.diff(y_major))), y_major[0]
    slope, intercept = np.polyfit(np.log10(X_LABELS), x_major, 1)
    resid = float(np.max(np.abs(np.polyval([slope, intercept], np.log10(X_LABELS)) - x_major)))
    if resid > 0.25:
        raise RuntimeError(f"x tick log fit residual {resid:.3f} pt -- axis may not be log")

    def to_data(x, y):
        return 10.0 ** ((x - intercept) / slope), 10.0 ** (Y_TOP_DECADE - (y - y_top) / y_per_dec)

    return to_data, {"x_per_decade": float(slope), "y_per_decade": y_per_dec, "x_resid": resid}


def extract(flavor):
    """(bottom, top, meta): the lower and upper boundaries of the low-mass BBN region."""
    page = fitz.open(str(require(BBN_HNL_PDFS[flavor])))[0]
    to_data, meta = calibrate(page)
    strokes, fills = [], []
    for d in page.get_drawings():
        if d["type"] == "s" and d["color"] and pdfio.color_matches(d["color"], NAVY) and not pdfio.is_dashed(d):
            strokes.append(pdfio.flatten(d["items"]))
        elif d["type"] in ("f", "fs") and d.get("fill") and pdfio.color_matches(d["fill"], NAVY):
            fills.append(pdfio.flatten(d["items"]))
    low = [s for s in strokes if max(to_data(x, y)[0] for x, y in s) < LOW_MASS_EDGE_GeV]
    region = [f for f in fills if max(to_data(x, y)[0] for x, y in f) < LOW_MASS_EDGE_GeV]
    if len(low) != 2 or len(region) != 1:
        raise RuntimeError(f"{flavor}: expected 2 boundary strokes and 1 fill below {LOW_MASS_EDGE_GeV} GeV, "
                           f"found {len(low)} and {len(region)}")
    fill = np.asarray(region[0])
    for stroke in low:
        gap = max(float(np.min(np.hypot(*(fill - p).T))) for p in stroke)
        if gap > 0.01:
            raise RuntimeError(f"{flavor}: a boundary stroke is {gap:.3f} pt off the fill outline")
    # Both strokes run left to right; the upper boundary ends higher (smaller PDF y).
    top_pts, bottom_pts = sorted(low, key=lambda s: s[-1][1])
    top = [to_data(x, y) for x, y in top_pts]
    bottom = [to_data(x, y) for x, y in bottom_pts]
    if any(np.diff([m for m, _ in top]) < 0) or any(np.diff([m for m, _ in bottom]) < 0):
        raise RuntimeError(f"{flavor}: a boundary is not single-valued in mass")
    meta.update(m_end=top[-1][0], u_top_end=top[-1][1],
                n_top=len(top), n_bottom=len(bottom))
    return bottom, top, meta


def crossing(curve, level):
    """Mass at which a descending curve first falls below ``level`` (log-log interpolation)."""
    for (m0, u0), (m1, u1) in zip(curve[:-1], curve[1:]):
        if u0 >= level > u1:
            t = (np.log10(level) - np.log10(u0)) / (np.log10(u1) - np.log10(u0))
            return 10.0 ** (np.log10(m0) + t * (np.log10(m1) - np.log10(m0)))
    return float("nan")


def write(flavor, bottom, top, meta):
    ylab, fig = PANELS[flavor]
    fname = BBN_HNL_PDFS[flavor].rsplit("/", 1)[1]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / f"BBN_ChenZhang2025_{flavor}.dat"
    lines = [
        f"BBN exclusion of a single Majorana HNL, {ylab} single-flavour mixing",
        f"Source: Chen and Zhang, 'BBN constraint on heavy neutrino production and decay',",
        f"  {ARXIV} (Phys. Rev. D 111 (2025) 123024), {fig}: the blue 'BBN Constraint' region.",
        f"Source figure PDF: sources/arXiv-2410.07343/{fname} (arXiv source bundle)",
        "Model: one Majorana heavy neutrino mixing with a single active flavour, whose",
        "  abundance (thermal freeze-out, or freeze-in from zero) and decays are both fixed",
        "  by (m_N, |U|^2) -- the paper's benchmark. Excluded at 95% CL (Y_p and D/H).",
        "  Used as published: no CL rescaling, no Dirac/Majorana factor.",
        "Method: vector extraction of the region's two solid navy boundary strokes,",
        "  RGB(0,0,0.502), each checked to lie on the outline of its 30%-opacity fill;",
        "  log-axis calibration from TICK GEOMETRY:",
        f"  x {meta['x_per_decade']:.4f} pt/decade (major ticks 1e-3 .. 10, 30 GeV; fit residual "
        f"{meta['x_resid']:.3f} pt),",
        f"  y {meta['y_per_decade']:.4f} pt/decade (a major tick every decade, 10^0 .. 10^-28).",
        f"Coverage: m_N = {top[0][0]:.3g} .. {meta['m_end']:.3g} GeV (top), {bottom[0][0]:.3g} .. "
        f"{bottom[-1][0]:.3g} GeV (bottom).",
        f"  The paper's region ends at m_N = {meta['m_end']:.3g} GeV with a dashed closure (top edge there",
        f"  {meta['u_top_end']:.3g}): its calculation stops, leaving 1-10 GeV as a gap. That closure",
        "  is not stored, nor are the >= 10 GeV region and, above 0.3 GeV, the CMB constraint,",
        "  both below |U|^2 = 1e-15.",
        f"  The top edge crosses |U|^2 = 1e-10 at m_N = {crossing(top, 1.0e-10):.3g} GeV.",
        "  The paper draws each boundary as a coarse polygon: the top edge has "
        f"{sum(0.3 <= m for m, _ in top)} vertices",
        f"  from 0.3 to {meta['m_end']:.3g} GeV, and is stored as drawn.",
        "top = upper boundary (the HNL decays too early, tau ~ 0.02 s); bottom = lower",
        "  boundary (too little abundance). The region between them is excluded.",
        f"Extracted: {DATE}",
        "mass_GeV  U2",
    ]
    with path.open("w") as fh:
        fh.writelines(f"# {line}\n" for line in lines)
        fh.write("# bottom\n" + "".join(f"{m:.6e}  {u:.6e}\n" for m, u in bottom))
        fh.write("\n# top\n" + "".join(f"{m:.6e}  {u:.6e}\n" for m, u in top))
    print(f"wrote {path}  (bottom {len(bottom)}, top {len(top)} vertices)")
    return path


def main() -> int:
    for flavor in PANELS:
        bottom, top, meta = extract(flavor)
        write(flavor, bottom, top, meta)
        print(f"   {flavor}: top edge {top[0][1]:.3g} at {top[0][0]:.3g} GeV .. {top[-1][1]:.3g} at "
              f"{top[-1][0]:.3g} GeV; crosses 1e-10 at {crossing(top, 1.0e-10):.3g} GeV")
    return 0


if __name__ == "__main__":
    sys.exit(main())
