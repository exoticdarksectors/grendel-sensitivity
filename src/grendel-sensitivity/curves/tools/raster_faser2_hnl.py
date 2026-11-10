"""FASER2 HNL projections (BC6, BC7, BC8) of Feng, Hewitt, Kling and La Rocco.

arXiv:2405.07330v2 (Phys. Rev. D 110, 035029), Fig. 9: the FASER2 (HL-LHC)
central line for |U_e|^2, |U_mu|^2 and |U_tau|^2, with the current FPF
geometry (650 m from the IP, 3 m x 1 m x 10 m), 3 ab^-1, Majorana HNLs from
pi, K, D, D_s, B, B_s, B_c and tau decays. It replaces the 2018 strawman
curves of arXiv:1811.12522 (HNLimits FASER_2018), which miss kaon production.

The figure is three bitmaps (2639 x 1085 px), so the route is raster:

* frame and tick marks from the dark pixels; y: major ticks at every decade
  (1e-3 at the top, 1e-8 at the bottom), fitted on the four inner ones and
  checked against the frame; x: one major tick (1 GeV), calibrated on the
  minor ticks 0.2-0.9 and 2-3 GeV fitted on log10(m), and checked against the
  frame's left edge (0.1 GeV);
* the light-blue central line (not its uncertainty band) isolated by colour;
  the "FASER2 (HL-LHC)" label shares the colour and is dropped as small (< 60 px)
  connected components;
* the line's centre points from short pixel runs, column- and row-wise (so
  steep parts are sampled too), chained by nearest neighbour from its left
  end and thinned to one point per 4 px.

The line is about 7 px wide, i.e. 0.04 decades in |U|^2: the digitisation is
good to a few per cent.

    python -m curves.tools.raster_faser2_hnl
"""
from __future__ import annotations

import sys

import fitz
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from ..plot.data import HNL_DATA
from . import pdfio
from .sources import FASER2_HNL_PDF, FASER2_HNL_PAGE, require

OUT_DIR = HNL_DATA / "raster"
DATE = "2026-10-25"
FLAVORS = ("Ue", "Umu", "Utau")    # the three bitmaps of Fig. 9, top to bottom
Y_TOP, Y_BOTTOM = -3.0, -8.0       # log10 |U|^2 at the frame's top and bottom edges


def bitmaps():
    doc = fitz.open(str(require(FASER2_HNL_PDF)))
    page = doc[FASER2_HNL_PAGE]
    placed = sorted((page.get_image_rects(xref)[0].y0, xref) for xref, *_ in page.get_images())
    if len(placed) != 3:
        raise RuntimeError(f"expected three bitmaps on the Fig. 9 page, found {len(placed)}")
    for _, xref in placed:
        pix = fitz.Pixmap(doc, xref)
        yield np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, pix.n)[:, :, :3].astype(int)


def groups(idx):
    out = []
    for i in idx:
        if out and i - out[-1][-1] <= 1:
            out[-1].append(i)
        else:
            out.append([i])
    return out


def calibration(img):
    dark = img.sum(axis=2) < 200
    h, w = dark.shape
    rows = np.where(dark.sum(axis=1) > 0.4 * w)[0]
    cols = np.where(dark.sum(axis=0) > 0.4 * h)[0]
    # Inner faces of the frame lines.
    top, bottom = groups(rows)[0][-1], groups(rows)[-1][0]
    left, right = groups(cols)[0][-1], groups(cols)[-1][0]
    # y ticks point into the frame from its left edge; majors are the longest.
    lengths = {}
    for j in range(top + 8, bottom - 8):
        run = dark[j, left + 1:left + 60]
        n = int(np.argmin(run)) if not run.all() else 59
        if n >= 6:
            lengths[j] = n
    ticks = [(float(np.mean(g)), max(lengths[j] for j in g)) for g in groups(sorted(lengths))]
    longest = max(n for _, n in ticks)
    ymaj = [y for y, n in ticks if n > 0.75 * longest]
    if len(ymaj) != 4:
        raise RuntimeError(f"expected four inner major y ticks, found {ymaj}")
    ypd = float(np.mean(np.diff(ymaj)))
    y0 = ymaj[0] - ypd    # 1e-3
    if abs(y0 - top) > 3 or abs(ymaj[-1] + ypd - bottom) > 3 or float(np.ptp(np.diff(ymaj))) > 2:
        raise RuntimeError(f"y ticks inconsistent with the frame: {ymaj}, frame {top}-{bottom}")
    # x ticks point up from the bottom edge.
    lengths = {}
    for i in range(left + 8, right - 8):
        run = dark[bottom - 1:bottom - 60:-1, i]
        n = int(np.argmin(run)) if not run.all() else 59
        if n >= 6:
            lengths[i] = n
    ticks = [(float(np.mean(g)), max(lengths[i] for i in g)) for g in groups(sorted(lengths))]
    longest = max(n for _, n in ticks)
    major = [x for x, n in ticks if n > 0.75 * longest]
    if len(major) != 1:
        raise RuntimeError(f"expected one major x tick (1 GeV), found {major}")
    # A tick may be hidden under a curve: each one is assigned the nearest k x 10^d of the provisional axis
    # through the frame's left edge (0.1 GeV) and the major tick, then all are fitted.
    pos = [x for x, _ in ticks]
    grid = np.log10([k * 10.0 ** d for d in (-1, 0) for k in range(1, 10)])
    rough = -1.0 + (np.asarray(pos) - left) / (major[0] - left)
    values = [float(10 ** grid[np.argmin(np.abs(grid - v))]) for v in rough]
    slope, icpt = np.polyfit(np.log10(values), pos, 1)
    resid = float(np.abs(np.polyval([slope, icpt], np.log10(values)) - pos).max())
    if resid > 1.5 or abs(np.polyval([slope, icpt], -1.0) - left) > 3:
        raise RuntimeError(f"x ticks off a log axis (residual {resid:.2f} px) or frame mismatch")

    def to_data(x, y):
        return 10 ** ((x - icpt) / slope), 10 ** (Y_TOP - (y - y0) / ypd)

    return to_data, (left, right, top, bottom), (slope, ypd, resid, values)


def line_points(img, frame):
    left, right, top, bottom = frame
    r, g, b = img[..., 0], img[..., 1], img[..., 2]
    mask = (r < 110) & (g > 185) & (g < 240) & (b > 235)
    mask[:top + 2] = mask[bottom - 1:] = False
    mask[:, :left + 2] = mask[:, right - 1:] = False
    lab, n = ndimage.label(mask, structure=np.ones((3, 3)))
    sizes = ndimage.find_objects(lab)
    keep = [k + 1 for k, s in enumerate(sizes) if max(s[0].stop - s[0].start, s[1].stop - s[1].start) > 60]
    mask = np.isin(lab, keep)
    pts = []
    for i in range(mask.shape[1]):
        for run in groups(np.where(mask[:, i])[0]):
            if len(run) <= 12:
                pts.append((i, float(np.mean(run))))
    for j in range(mask.shape[0]):
        for run in groups(np.where(mask[j])[0]):
            if len(run) <= 12:
                pts.append((float(np.mean(run)), j))
    return np.asarray(pts, float)


def chain(points, step=4.0, reach=80.0):
    # Other limits drawn over the line cut it into pieces up to ~50 px apart; the next point is always the
    # nearest free one, so a gap is bridged to the continuation, not to the other branch.
    tree = cKDTree(points)
    free = np.ones(len(points), bool)
    cur = int(np.argmin(points[:, 0] + 1e-3 * -points[:, 1]))   # leftmost, lowest on screen if tied
    out = [points[cur]]
    while True:
        free[tree.query_ball_point(points[cur], step)] = False
        near = [k for k in tree.query_ball_point(points[cur], reach) if free[k]]
        if not near:
            break
        cur = min(near, key=lambda k: np.hypot(*(points[k] - points[cur])))
        out.append(points[cur])
    if free.mean() > 0.02:
        raise RuntimeError(f"{free.mean():.1%} of the line points were not reached by the chain")
    return out


def main() -> int:
    for flavor, img in zip(FLAVORS, bitmaps()):
        to_data, frame, (slope, ypd, resid, values) = calibration(img)
        ordered = chain(line_points(img, frame))
        lo, hi = to_data(frame[0], frame[2]), to_data(frame[1], frame[3])
        pdfio.write_chains(OUT_DIR / f"FASER2_2405_{flavor}.dat", [[to_data(x, y) for x, y in ordered]],
                           [f"FASER2 (HL-LHC) HNL |{flavor[1:]}|^2 sensitivity, Majorana, 3 ab^-1, FPF geometry",
                            "Digitised from arXiv:2405.07330v2 (Phys. Rev. D 110, 035029), Fig. 9, central line",
                            "(not the uncertainty band); raster route, see curves/tools/raster_faser2_hnl.py.",
                            f"Date: {DATE}",
                            f"Calibration: x {slope:.2f} px/decade fitted on the ticks {values[0]:.1f}-{values[-1]:.0f} GeV",
                            f"(residual {resid:.2f} px), y {ypd:.2f} px/decade (major ticks at every decade).",
                            f"Frame m_N in [{lo[0]:.3g}, {hi[0]:.3g}] GeV, |U|^2 in [{hi[1]:.3g}, {lo[1]:.3g}];",
                            "the line leaves the frame at its top (1e-3)."],
                           columns='Format: chains separated by blank line; each row "mass_GeV  U2"')
    return 0


if __name__ == "__main__":
    sys.exit(main())
