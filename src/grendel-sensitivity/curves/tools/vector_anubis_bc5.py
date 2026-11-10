"""ANUBIS BC5 (BR(h -> SS) = 0.01) from the collaboration's own paper.

arXiv:2510.26932v2 (doi:10.1103/kd5j-d4fv), Fig. 14a,
``figures/Combined_CLsOnly/sinTheta/sin_theta_br0p01_SigEff0p5_ANUBISonly_ceilingOnly_withBand_CombOnly.pdf``:
ceiling configuration, 1+2DVs, 95% CL exclusion, eps_sig = 0.5, 3 ab^-1,
m_S = 10-60 GeV. Two closed contours of the same edge colour: the
conservative-background one (outline over a hatch pattern) and the
background-free one (filled). The conservative contour replaces the PBC
2025 ANUBIS curve (``anubis_pbc_bc5.dat``, kept for reference) in the
paper panel: that is the collaboration's ESPP input, superseded by this
paper (its earlier conservative projections used 5 sigma discovery, this
one 95% CL), and the Higgs-plane figure of the paper (Fig. 3) already
shows the conservative 1DV projection of the same paper. The
background-free contour is written for reference only.

Calibration from the frame: x linear, 97.608 -> 696.510 pt = 10 -> 60 GeV;
y log, 511.864 -> 19.404 pt = 1e-10 -> 1e-3 in sin(theta) (major ticks
every 70.351 pt). Values are squared to sin^2(theta).

    python -m curves.tools.vector_anubis_bc5
"""
from __future__ import annotations

import sys

import fitz

from ..plot.data import BC5_DATA
from . import pdfio
from .sources import ANUBIS_BC5_PDF, require
from .vector_bc4 import write

PROJ_DIR = BC5_DATA / "projections"
DATE = "2026-10-25"
EDGE = (0.0, 0.674, 0.9)
X0, X1, M0, M1 = 97.608, 696.510, 10.0, 60.0
Y_BOTTOM, Y_TOP, DECADES = 511.864, 19.404, 7.0
CONTOURS = {
    # file stem, path type, scenario
    "anubis_conservative_bc5": ("s", "conservative background (hatched in the source); drawn in the paper panel"),
    "anubis_bkgfree_bc5": ("fs", "background-free (filled in the source); reference only, not drawn"),
}


def anubis_to_data(x, y):
    log_sinth = -10.0 + (Y_BOTTOM - y) / (Y_BOTTOM - Y_TOP) * DECADES
    return M0 + (x - X0) / (X1 - X0) * (M1 - M0), 10 ** (2.0 * log_sinth)


def main() -> int:
    page = fitz.open(str(require(ANUBIS_BC5_PDF)))[0]
    for stem, (kind, scenario) in CONTOURS.items():
        # The hatch-pattern tiles and legend patches share the edge colour; only
        # the contours span the frame from 10 to 60 GeV.
        paths = [d for d in page.get_drawings() if d["type"] == kind and d.get("color")
                 and pdfio.color_matches(d["color"], EDGE, 5e-3)
                 and abs(d["rect"].x0 - X0) < 1 and abs(d["rect"].x1 - X1) < 1]
        if len(paths) != 1:
            raise RuntimeError(f"{stem}: expected one contour, found {len(paths)}")
        items = paths[0]["items"]
        chain = [(it[1].x, it[1].y) for it in items] + [(items[-1][2].x, items[-1][2].y)]
        write(PROJ_DIR / f"{stem}.dat", [chain], anubis_to_data,
              ["ANUBIS BC5 dark-scalar sensitivity, BR(h -> SS) = 0.01",
               "Extracted from arXiv:2510.26932v2 / Fig. 14a, "
               "figures/Combined_CLsOnly/sinTheta/sin_theta_br0p01_SigEff0p5_ANUBISonly_ceilingOnly_withBand_CombOnly.pdf",
               f"Date: {DATE}",
               f"Ceiling, 1+2DVs, 95% CL exclusion, eps_sig = 0.5, 3 ab^-1; {scenario}.",
               "Closed contour with vertices at 10, 15, 20, 25, 30, 40, 50, 60 GeV.",
               "Source y axis is sin(theta); values below are squared to sin^2(theta)."])
    return 0


if __name__ == "__main__":
    sys.exit(main())
