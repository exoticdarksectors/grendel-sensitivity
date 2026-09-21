"""BC7 (HNL, muon dominance) from the PBC 2025 figure, vector route.

``pbc_bc7_gray_final.pdf`` of arXiv:2505.00947v2 (label
``fig:pbc_bc7_gray_final``). Two extractions:

* ``islands``: for each contour-coloured stroke (ANUBIS, CODEX-b, FLArE,
  FASER2, SHiP) the line items are stitched into a polyline, transformed
  into (m_N, U^2) with the axis-frame calibration, and sampled on a mass
  grid: every intersection with the vertical line m = m_sample, smallest
  and largest U^2, gives the island's ``u2_min``/``u2_max`` at that mass.
* ``excluded``: the gray "currently excluded" fill, a single drawing with
  ~200 path items holding several disjoint subpaths (one per experimental
  limit merged into the union). Each subpath is written as its own polygon.

    python -m curves.tools.vector_bc7 [islands|excluded|all]
"""
from __future__ import annotations

import sys

import fitz
import numpy as np

from ..plot.data import HNL_DATA
from . import pdfio
from .sources import BC7_PDF, require

OUT_DIR = HNL_DATA / "vector"
SOURCE = "Figures/5-PhysicsReach/BSM_benchmarks/pbc_bc7_gray_final.pdf"
FIG_LABEL = "fig:pbc_bc7_gray_final"
DATE = "2026-06-06"

# Axis-frame calibration from the four black stroke lines.
X_LEFT, X_RIGHT = 79.69999694824219, 526.0999755859375
Y_BOT, Y_TOP = 359.79840087890625, 27.158416748046875
LOG_M_LEFT, LOG_M_RIGHT = -1.0, 2.0
LOG_U_BOT, LOG_U_TOP = -12.0, -2.0

CURVES = {
    "ANUBIS": (0.29411765933036804, 0.0, 0.5098039507865906),
    "CODEX-b": (0.0, 0.545098066329956, 0.545098066329956),
    "FLArE": (1.0, 0.0, 0.0),
    "FASER2": (0.0, 0.3921568691730499, 0.0),
    "SHiP": (0.8549019694328308, 0.43921568989753723, 0.8392156958580017),
}
GRAY = (0.5018997192382812, 0.5019760727882385, 0.5018844604492188)
SENTINEL = 1.0e-1


def pdf_to_data(x, y):
    log_m = LOG_M_LEFT + (x - X_LEFT) / (X_RIGHT - X_LEFT) * (LOG_M_RIGHT - LOG_M_LEFT)
    log_u2 = LOG_U_BOT + (Y_BOT - y) / (Y_BOT - Y_TOP) * (LOG_U_TOP - LOG_U_BOT)
    return 10 ** log_m, 10 ** log_u2


def build_mass_grid():
    # Dense log-uniform spine + finer steps across the kaon/D/B kinematic features.
    base = np.logspace(np.log10(0.1), np.log10(30), 300)
    dense = np.concatenate([np.arange(0.25, 0.45 + 1e-9, 0.002), np.arange(1.60, 2.00 + 1e-9, 0.002),
                            np.arange(4.80, 5.40 + 1e-9, 0.005)])
    return np.sort(np.unique(np.round(np.concatenate([base, dense]), 5)))


MASS_GRID = build_mass_grid()


def per_curve_grid(poly_data):
    """The global grid plus every native polyline x: the exact masses where
    the figure's polyline has a kink, so sampling there keeps every feature."""
    native = sorted({round(m, 6) for m, _ in poly_data})
    return np.sort(np.unique(np.concatenate([MASS_GRID, native])))


def vertical_crossings(poly_data, m_target, eps=1e-9):
    """All U^2 values where the polyline crosses the vertical line m = m_target."""
    crossings = []
    lm_target = np.log10(m_target)
    for i in range(len(poly_data) - 1):
        m1, u1 = poly_data[i]
        m2, u2 = poly_data[i + 1]
        lm1, lm2 = np.log10(m1), np.log10(m2)
        if (lm1 - lm_target) * (lm2 - lm_target) <= 0:
            if abs(lm2 - lm1) < eps:
                crossings.append(u1)
                crossings.append(u2)
            else:
                t = (lm_target - lm1) / (lm2 - lm1)
                lu = np.log10(u1) + t * (np.log10(u2) - np.log10(u1))
                crossings.append(10 ** lu)
    return crossings


def island_rows(name, poly):
    m_min_curve = min(p[0] for p in poly)
    m_max_curve = max(p[0] for p in poly)
    rows = []
    for m in per_curve_grid(poly):
        if m < m_min_curve - 1e-6 or m > m_max_curve + 1e-6:
            continue
        cr = vertical_crossings(poly, m)
        if not cr:
            continue
        u_lo, u_hi = min(cr), max(cr)
        if name == "SHiP" and u_hi > 1.0e-2:
            u_hi = SENTINEL
        rows.append((float(m), float(u_lo), float(u_hi)))
    return rows


def write_island(name, fname, rows, extra=None):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / fname
    with path.open("w") as fh:
        fh.write(f"# {name} HNL |Umu|^2 sensitivity (BC7 benchmark)\n")
        fh.write(f"# Digitized from arXiv:2505.00947v2, {FIG_LABEL}, {DATE}\n")
        fh.write("# Method: vector extraction from PDF path operators (PyMuPDF)\n")
        fh.write(f"# Source file: {SOURCE}\n")
        for ln in extra or ():
            fh.write(f"# {ln}\n")
        fh.write("# mass_GeV  u2_min  u2_max\n")
        for m, lo, hi in rows:
            fh.write(f"{m:.4f}  {lo:.3e}  {hi:.3e}\n")
    print(f"wrote {path} ({len(rows)} rows)")


def islands(page):
    draws = page.get_drawings()
    polys = {}
    for name, target in CURVES.items():
        matches = [d for d in draws if d.get("type") == "s" and d.get("color") and pdfio.near(d["color"], target)]
        matches.sort(key=lambda d: -len(d.get("items", [])))
        pts = pdfio.stitch(matches[0]["items"])
        # The arrow-leader Bezier is a separate single-item drawing, never the longest match.
        polys[name] = [pdf_to_data(x, y) for x, y in pts]
        print(f"{name}: {len(pts)} vertices, m in [{min(p[0] for p in polys[name]):.3g}, "
              f"{max(p[0] for p in polys[name]):.3g}], U2 in [{min(p[1] for p in polys[name]):.2e}, "
              f"{max(p[1] for p in polys[name]):.2e}]")
    rows = {name: island_rows(name, poly) for name, poly in polys.items()}
    overlap = ["FLArE and FASER2 contours overlap closely; vector paths give two",
               "independent extractions and reveal whether/where they differ."]
    write_island("ANUBIS", "ANUBIS_Umu.dat", rows["ANUBIS"])
    write_island("CODEX-b", "CODEX-b_Umu.dat", rows["CODEX-b"])
    write_island("FLArE", "FLArE_Umu.dat", rows["FLArE"], extra=overlap)
    write_island("FASER2", "FASER2_Umu.dat", rows["FASER2"],
                 extra=["FASER2 and FLArE contours overlap closely; vector paths give two",
                        "independent extractions and reveal whether/where they differ."])
    write_island("SHiP", "SHiP_Umu.dat", rows["SHiP"],
                 extra=["Where the SHiP contour upper crossing exceeds the plotted ceiling",
                        "(U²>1e-2), the upper edge is clipped to sentinel value 1.0e-1."])


def excluded(page):
    main = next(d for d in page.get_drawings() if d.get("type") in ("f", "fs") and d.get("fill")
                and pdfio.near(d["fill"], GRAY) and len(d.get("items", [])) > 100)
    subs = pdfio.subpaths(main["items"])
    print(f"{len(subs)} subpaths, sizes: {[len(s) for s in subs[:10]]} ...")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / "currently_excluded.dat"
    total = 0
    with path.open("w") as fh:
        fh.write('# "Currently excluded" gray region from BC7 fig\n')
        fh.write(f"# Extracted from arXiv:2505.00947v2 / {SOURCE}\n")
        fh.write(f"# Date: {DATE}\n")
        fh.write('# Format: subpaths separated by blank line; each row is "mass_GeV  U2"\n')
        fh.write("# Source experiments per caption: NuTeV, BEBC, T2K, NA62, E949,\n")
        fh.write("# MicroBooNE, CMS, ATLAS, CHARM-II, DELPHI, KOTO, LHCb, BelleII, plus\n")
        fh.write("# reinterpretations of PS191, CHARM, LSND.\n")
        fh.write("# mass_GeV  U2\n")
        for k, sub in enumerate(subs):
            if k > 0:
                fh.write("\n")
            fh.write(f"# subpath {k}: {len(sub)} vertices\n")
            for x, y in sub:
                m, u2 = pdf_to_data(x, y)
                if m < 1e-3 or m > 1e3 or u2 < 1e-14:
                    continue
                fh.write(f"{m:.5e}  {u2:.5e}\n")
                total += 1
    print(f"wrote {path} ({total} vertices total)")


def main(argv=None) -> int:
    what = (argv if argv is not None else sys.argv[1:]) or ["all"]
    page = fitz.open(str(require(BC7_PDF)))[0]
    if what[0] in ("islands", "all"):
        islands(page)
    if what[0] in ("excluded", "all"):
        excluded(page)
    return 0


if __name__ == "__main__":
    sys.exit(main())
