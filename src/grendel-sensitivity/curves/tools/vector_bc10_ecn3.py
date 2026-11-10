"""SHiP BC10 curve of the PBC ECN3 report figure (a reference, not drawn since
``senscalc_bc10_ship``; it used to fill in below the eta).

``Figures/9-PhysicsPotential/BC10.pdf`` of arXiv:2310.17726 (g_Y = 2v/f_a,
m 0.03-0.975 GeV, g_Y 1e-6..1e-3): SHiP is the red dashed curve. Where it
overlaps the arXiv:2501.04525 curve (0.62-0.82 GeV) the two agree within
3 % in g_Y, so it carries the same production (the B -> K^(i) tower); it
continues below the eta, where arXiv:2501.04525 has none. The axes are
calibrated on the 1e-1 GeV and the 1e-6..1e-3 major ticks, with the
x scale from the 0.2 GeV minor tick (every minor tick within 0.05 pt).

Output: ``curves/data/bc10/projections/ship_ecn3.dat``, rows
``mass_GeV  invf_GeV^-1`` (1/f BNT = g_Y / v_h).

    python -m curves.tools.vector_bc10_ecn3
"""
from __future__ import annotations

import sys

import fitz
import numpy as np

from ..plot.data import BC10_DATA
from . import pdfio
from .sources import ECN3_BC10_PDF, require

V_H = 246.0
DATE = "2026-10-07"
X_M1, X_PER_DEC = 222.59, (311.47 - 222.59) / np.log10(2.0)   # 0.1 GeV major tick, 0.2 GeV minor tick
Y_M6, Y_PER_DEC = 345.39, (345.39 - 12.75) / 3.0               # 1e-6 and 1e-3 major ticks


def to_data(x, y):
    m = 10 ** (-1.0 + (x - X_M1) / X_PER_DEC)
    gy = 10 ** (-6.0 + (Y_M6 - y) / Y_PER_DEC)
    return m, gy / V_H


def main() -> int:
    page = fitz.open(str(require(ECN3_BC10_PDF)))[0]
    for m, x in ((0.03, 68.2), (0.3, 363.46), (0.5, 428.96), (0.9, 504.34)):   # minor ticks
        if abs(X_M1 + np.log10(m / 0.1) * X_PER_DEC - x) > 0.05:
            raise RuntimeError(f"x calibration off at the {m} GeV tick")
    chains = pdfio.collect_strokes(page, (1.0, 0.0, 0.0), t_samples=pdfio.linspace_samples(8), tol=0.01)
    pdfio.write_chains(
        BC10_DATA / "projections" / "ship_ecn3.dat", [[to_data(x, y) for x, y in ch] for ch in chains],
        ["SHiP fermiophilic-ALP (BC10) sensitivity, PBC ECN3 report (arXiv:2310.17726, Fig. BC10)",
         f"Date: {DATE}",
         "Source y axis g_Y = 2 v/f_a; stored as 1/f (BNT) = g_Y/246. The frame ends at 0.975 GeV;",
         "the path beyond it was not visually published. Within 3 % of arXiv:2501.04525 at 0.62-0.82 GeV."],
        columns='Format: chains separated by blank line; each row "mass_GeV  invf_GeV^-1"')
    return 0


if __name__ == "__main__":
    sys.exit(main())
