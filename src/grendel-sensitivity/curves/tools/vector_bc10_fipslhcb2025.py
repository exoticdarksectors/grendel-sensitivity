"""BC10 layers of the FIPs at LHCb 2025 summary figure.

``ALPs_fermions_FIPsLHCb25_all.pdf`` of arXiv:2510.05257v2 (g_Y = 2v/f vs
m_ALP, m 1e-2..20 GeV, g_Y 1e-6..1e-1) is the most recent BC10 compilation.
Its gray regions are the FIPs 2022 polygons (``vector_bc10_fips2022``,
identical below g_Y = 0.1); new are the NA62 regions and the projections:

    fill (0.2,0.2,0.4) opacity 0.4   NA62: beam-dump 2025 (arXiv:2502.04241),
                                     K+ -> pi+ X and pi0 -> invisible
    blue long dash                   CODEX-b (arXiv:2505.05952)
    magenta dotted                   FASER2 (arXiv:1811.12522)

Its SHiP curve (red long dash, the ECN3 LoI) is not used: it lies 1.41x
above the SHiP curve of the PBC ECN3 figure in g_Y, the B -> K^(i) tower
it lacks (``vector_bc10_ecn3``). The axes are calibrated on the major ticks.

Output: ``curves/data/bc10/vector/fipslhcb2025/na62_*.dat`` and
``curves/data/bc10/projections/{codexb,faser2}_fipslhcb25.dat``, rows
``mass_GeV  invf_GeV^-1`` (1/f BNT = g_Y / v_h).

    python -m curves.tools.vector_bc10_fipslhcb2025
"""
from __future__ import annotations

import sys

import fitz
import numpy as np

from ..plot.data import BC10_DATA
from . import pdfio
from .sources import FIPSLHCB2025_ALP_PDF, require

V_H = 246.0
DATE = "2026-10-07"
SAMPLES = pdfio.linspace_samples(6)
SOURCE = "Extracted from arXiv:2510.05257v2 (FIPs at LHCb 2025), ALPs_fermions_FIPsLHCb25_all.pdf"
COLUMNS = 'Format: chains separated by blank line; each row "mass_GeV  invf_GeV^-1"'


class Frame:
    """Log axes from the major ticks on the bottom and left frame edges:
    x majors 1e-2, 1e-1, 1, 10; y majors 1e-1 .. 1e-6."""

    def __init__(self, page):
        black = [d for d in page.get_drawings() if d["type"] == "s" and d["color"]
                 and pdfio.color_matches(d["color"], (0.0, 0.0, 0.0)) and not pdfio.is_dashed(d)]
        frame = max(black, key=lambda d: d["rect"].width * d["rect"].height)["rect"]
        xt, yt = [], []
        for d in black:
            for it in d["items"]:
                if it[0] != "l":
                    continue
                a, b = it[1], it[2]
                if abs(a.x - b.x) < 0.05 and min(abs(a.y - frame.y1), abs(b.y - frame.y1)) < 0.6:
                    xt.append((abs(a.y - b.y), a.x))
                if abs(a.y - b.y) < 0.05 and min(abs(a.x - frame.x0), abs(b.x - frame.x0)) < 0.6:
                    yt.append((abs(a.x - b.x), a.y))
        x_major = sorted({round(x, 2) for length, x in xt if 8.0 < length < 20.0})
        y_major = sorted({round(y, 2) for length, y in yt if 10.0 < length < 20.0})
        if len(x_major) != 4 or len(y_major) != 6:
            raise RuntimeError(f"unexpected major ticks: x {x_major}, y {y_major}")
        self.x0, self.xd = x_major[0], float(np.mean(np.diff(x_major)))   # 1e-2, pt per decade
        self.y0, self.yd = y_major[0], float(np.mean(np.diff(y_major)))   # 1e-1, pt per decade
        print(f"x: 1e-2 at {self.x0}, {self.xd:.3f} pt/decade; y: 1e-1 at {self.y0}, {self.yd:.3f} pt/decade")

    def to_data(self, x, y):
        m = 10 ** (-2.0 + (x - self.x0) / self.xd)
        gy = 10 ** (-1.0 - (y - self.y0) / self.yd)
        return m, gy / V_H


def strokes(page, color, dashes_prefix):
    chains = [pdfio.flatten(d["items"], SAMPLES) for d in page.get_drawings()
              if d["type"] == "s" and d["color"] and pdfio.color_matches(d["color"], color, 0.01)
              and str(d.get("dashes") or "").startswith(dashes_prefix)]
    return pdfio.merge_chains([c for c in chains if len(c) >= 2])


def write(path, frame, chains, lines):
    pdfio.write_chains(path, [[frame.to_data(x, y) for x, y in ch] for ch in chains],
                       [*lines, SOURCE, f"Date: {DATE}", "Source y axis g_Y = 2 v_h/f; stored as 1/f (BNT) = g_Y/246."],
                       columns=COLUMNS)


def main() -> int:
    page = fitz.open(str(require(FIPSLHCB2025_ALP_PDF)))[0]
    frame = Frame(page)
    na62 = {}
    for ch in pdfio.collect_fills(page, (0.2, 0.2, 0.4), t_samples=SAMPLES, opacity=0.4, rects=True):
        m = [frame.to_data(x, y)[0] for x, y in ch]
        name = "na62_pi0_inv" if min(m) > 0.1 else ("na62_dump_2025" if max(m) > 0.4 else "na62_kpix")
        na62[name] = ch
    descriptions = {"na62_dump_2025": "NA62 beam-dump search, 2021 data (arXiv:2502.04241) -- existing exclusion",
                    "na62_kpix": "NA62 K+ -> pi+ X (invisible) -- existing exclusion",
                    "na62_pi0_inv": "NA62 pi0 -> invisible window -- existing exclusion"}
    for name, ch in na62.items():
        write(BC10_DATA / "vector" / "fipslhcb2025" / f"{name}.dat", frame, [ch], [descriptions[name]])
    proj = BC10_DATA / "projections"
    write(proj / "codexb_fipslhcb25.dat", frame, strokes(page, (0.0, 0.0, 1.0), "[ 19"),
          ["CODEX-b baseline (arXiv:2505.05952, 300 fb^-1); upper edge cut by the frame top g_Y = 0.1.",
           "Equal to that paper's alp_f.pdf (c_f at Lambda = 1 TeV, g_Y = v_h c_f/Lambda) to <1% at 0.5-3 GeV."])
    write(proj / "faser2_fipslhcb25.dat", frame, strokes(page, (1.0, 0.0, 1.0), "[ 2"),
          ["FASER2 (arXiv:1811.12522)"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
