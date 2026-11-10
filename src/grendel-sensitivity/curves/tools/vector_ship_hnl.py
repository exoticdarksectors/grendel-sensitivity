"""SHiP HNL projections at 6e20 protons on target (BC6, BC7, BC8).

The HNLimits SHiP curves (arXiv:1811.00930) are for 2e20 protons on target,
a third of the exposure of every other SHiP curve in the paper. The 6e20
curves:

* BC6 (|U_e|^2) and BC7 (|U_mu|^2): ESPP 2026 Physics Briefing Book,
  arXiv:2511.03883v1, Fig. 8.19 (page index 145): bottom panel |U_e|^2, top
  panel |U_mu|^2, royal-blue "SHiP, 15y, 6x10^20 PoT" curve. The BC7 curve
  equals the PBC 2025 one (``vector/SHiP_Umu.dat``). Both axes carry a major
  tick at every decade (m_N 0.1 - 1e4 GeV, |U|^2 1e-1 - 1e-12).
* BC8 (|U_tau|^2): the only 6e20 curve is that of the PBC ECN3 report,
  arXiv:2310.17726v1, ``Figures/9-PhysicsPotential/BC8.pdf``, red dashed
  SHiP curve; its upper edge runs out of the frame at |U_tau|^2 = 1e-3 and
  its lower edge starts at the frame's 0.2 GeV. The x axis has a single
  major tick (1 GeV): it is calibrated from the minor ticks 0.2-0.9 and
  2-5 GeV, fitted on log10(m) with a residual below 0.1 pt.

The 6e20 curves drop B_c production (their mass reach, 5.1 GeV for e and mu,
3.5 GeV for tau, is that of the no-B_c curves of arXiv:1811.00930 at three
times the exposure), so they reach smaller couplings but lower masses.

    python -m curves.tools.vector_ship_hnl
"""
from __future__ import annotations

import sys

import fitz
import numpy as np

from ..plot.data import HNL_DATA
from . import pdfio
from .sources import BB_HNL_PAGE, BB_PDF, ECN3_BC8_PDF, require
from .vector_bc4 import BB_SHIP_BLUE, ECN3_SHIP_RED, SAMPLES, _black_lines, major_y_calibration, \
    panel_tick_calibration

OUT_DIR = HNL_DATA / "vector"
DATE = "2026-10-25"
COLUMNS = 'Format: chains separated by blank line; each row "mass_GeV  U2"'


def panel_frames(page):
    """(left, right, top, bottom) of each plot frame on a page, top to bottom."""
    hs = sorted({(round(a.y, 2), round(min(a.x, b.x), 2), round(max(a.x, b.x), 2))
                 for a, b in _black_lines(page) if abs(a.y - b.y) < 0.05 and abs(a.x - b.x) > 50})
    frames = []
    while hs:
        top = hs.pop(0)
        bottom = next(h for h in hs if abs(h[1] - top[1]) < 0.1 and abs(h[2] - top[2]) < 0.1)
        hs.remove(bottom)
        frames.append((top[1], top[2], top[0], bottom[0]))
    return frames


def write(path, chains, to_data, header, clip_x, clip_y):
    pdfio.write_chains(path, pdfio.clip_chains(chains, to_data, clip_x, clip_y), header, columns=COLUMNS)


def briefing_book():
    page = fitz.open(str(require(BB_PDF)))[BB_HNL_PAGE]
    frames = panel_frames(page)
    if len(frames) != 2:
        raise RuntimeError(f"expected two panels on the Fig. 8.19 page, found {len(frames)}")
    for flavor, frame, panel in (("Umu", frames[0], "top"), ("Ue", frames[1], "bottom")):
        left, right, top, bottom = frame
        xmaj, xpd, ymaj, ypd = panel_tick_calibration(page, frame)
        if len(xmaj) != 6 or len(ymaj) != 12 or abs(xmaj[0] - left) > 0.1 or abs(ymaj[0] - top) > 0.1:
            raise RuntimeError(f"unexpected tick pattern: x {xmaj}, y {ymaj}")

        def to_data(x, y, xmaj=xmaj, xpd=xpd, ymaj=ymaj, ypd=ypd):
            return 10 ** (-1.0 + (x - xmaj[0]) / xpd), 10 ** (-1.0 - (y - ymaj[0]) / ypd)

        (m_lo, u_top), (m_hi, u_bottom) = to_data(left, top), to_data(right, bottom)
        chains = [sub for d in page.get_drawings()
                  if d["type"] == "s" and d["color"] and pdfio.color_matches(d["color"], BB_SHIP_BLUE)
                  and top - 1 < d["rect"].y0 and d["rect"].y1 < bottom + 1 and d["rect"].x1 < right + 1
                  for sub in pdfio.subpaths(d["items"], SAMPLES)]
        write(OUT_DIR / f"SHiP_bb_{flavor}.dat", pdfio.merge_chains(chains), to_data,
              [f"SHiP HNL |{flavor[1:]}|^2 sensitivity, Majorana, 15 y, 6e20 protons on target",
               f"Extracted from arXiv:2511.03883v1 (ESPP 2026 Physics Briefing Book), Fig. 8.19 {panel} panel,",
               "royal-blue SHiP curve (cite this).",
               f"Date: {DATE}",
               f"Calibration from tick geometry: x {xpd:.4f} pt/decade, y {ypd:.4f} pt/decade (major ticks at",
               "every decade, every minor tick on a log10(k) position).",
               f"Clipped to the frame m_N in [{m_lo:.3g}, {m_hi:.3g}] GeV, |U|^2 in [{u_bottom:.3g}, {u_top:.3g}]."],
              (m_lo, m_hi), (u_bottom, u_top))


def ecn3_bc8():
    page = fitz.open(str(require(ECN3_BC8_PDF)))[0]
    xs = [a.x for a, b in _black_lines(page) if abs(a.x - b.x) < 0.05 and abs(a.y - b.y) > 50]
    ys = [a.y for a, b in _black_lines(page) if abs(a.y - b.y) < 0.05 and abs(a.x - b.x) > 50]
    left, right, top, bottom = min(xs), max(xs), min(ys), max(ys)
    # The ticks point outwards: they start on the frame's bottom edge.
    marks = sorted({(round((a.x + b.x) / 2, 3), round(abs(a.y - b.y), 2)) for a, b in _black_lines(page)
                    if abs(a.x - b.x) < 0.05 and 0.5 < abs(a.y - b.y) < 25 and abs(min(a.y, b.y) - bottom) < 0.05
                    and left - 1 < a.x < right + 1})
    longest = max(w for _, w in marks)
    major = [x for x, w in marks if abs(w - longest) < 0.2]
    if len(major) != 1:
        raise RuntimeError(f"expected one major x tick (1 GeV), found {major}")
    below = [x for x, w in marks if x < major[0]]
    above = [x for x, w in marks if x > major[0]]
    values = [0.1 * k for k in range(10 - len(below), 10)] + [1.0] + [float(k) for k in range(2, 2 + len(above))]
    pos = below + major + above
    slope, icpt = np.polyfit(np.log10(values), pos, 1)
    resid = np.abs(np.polyval([slope, icpt], np.log10(values)) - pos)
    if resid.max() > 0.1:
        raise RuntimeError(f"x ticks off a log axis by up to {resid.max():.3f} pt: {list(zip(values, pos))}")
    cy, ypd = major_y_calibration(page, decades_between_major=1.0)
    top_label = -3.0  # the top major tick is 1e-3

    def to_data(x, y):
        return 10 ** ((x - icpt) / slope), 10 ** (top_label - (y - cy) / ypd)

    (m_lo, u_top), (m_hi, u_bottom) = to_data(left, top), to_data(right, bottom)
    chains = [sub for d in page.get_drawings()
              if d["type"] == "s" and d["color"] and pdfio.color_matches(d["color"], ECN3_SHIP_RED)
              and pdfio.is_dashed(d) for sub in pdfio.subpaths(d["items"], SAMPLES)]
    write(OUT_DIR / "SHiP_ecn3_Utau.dat", pdfio.merge_chains(chains), to_data,
          ["SHiP HNL |U_tau|^2 sensitivity, 6e20 protons on target",
           "Extracted from arXiv:2310.17726v1 (PBC, Post-LS3 Experimental Options in ECN3),",
           "Figures/9-PhysicsPotential/BC8.pdf, red dashed SHiP curve (cite this).",
           f"Date: {DATE}",
           f"Calibration from tick geometry: x {slope:.4f} pt/decade fitted on the ticks {values[0]:.1f}-{values[-1]:.0f} GeV",
           f"(residual {resid.max():.3f} pt), y {ypd:.4f} pt/decade (major ticks at every decade, 1e-3 at the top).",
           f"Clipped to the frame m_N in [{m_lo:.3g}, {m_hi:.3g}] GeV, |U|^2 in [{u_bottom:.3g}, {u_top:.3g}]:",
           "the upper edge leaves the frame at 1e-3 and the lower edge starts at its left edge."],
          (m_lo, m_hi), (u_bottom, u_top))


def main() -> int:
    briefing_book()
    ecn3_bc8()
    return 0


if __name__ == "__main__":
    sys.exit(main())
