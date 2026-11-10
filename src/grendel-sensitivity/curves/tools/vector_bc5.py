"""BC5 (dark scalar, mixing plus quartic, BR(h -> SS) = 0.01) reference curves.

All from the PBC 2025 BC5 figure (arXiv:2505.00947v2,
``pbc_bc5_gray_final.pdf``; its arXiv source bundle carries no data files,
so the curves are read from the figure's path operators). The existing
exclusion of that figure is mixing-only and is shared with BC4
(``curves/data/bc4/vector/currently_excluded.dat``, ``vector_bc4``).

Projections kept (dedicated LLP experiments only): ANUBIS, CODEX-b,
MATHUSLA40, SHiP. FASER2, FLArE, MAPP, MAPP2 and the LHCb downstream
recast are not extracted. Each curve is the width-2 stroke of its colour;
the width-1 stroke of the same colour is the label leader line and is
dropped. Per the caption the line style is the maturity of the
background estimate: solid = extrapolation of existing data (ANUBIS),
dashed = full simulation (MATHUSLA40, SHiP), dotted = toy simulation or
negligible background (CODEX-b). Axis calibration and the sin(theta) ->
sin^2(theta) conversion are ``vector_bc4.pbc_to_data``.

    python -m curves.tools.vector_bc5
"""
from __future__ import annotations

import sys

import fitz

from ..plot.data import BC5_DATA
from . import pdfio
from .sources import BC5_PDF, require
from .vector_bc4 import BC5_SOURCE, SAMPLES, pbc_to_data, write

PROJ_DIR = BC5_DATA / "projections"
DATE = "2026-09-26"
CURVES = {
    # file stem, stroke colour, caption line style, extra notes
    "anubis": ((0.294, 0.0, 0.51), "solid (background from extrapolation of existing data)",
               ["The published contour is a 16-segment polyline from 5 to 62.5 GeV."]),
    "codexb": ((0.0, 0.545, 0.545), "dotted (toy simulation or negligible background)", []),
    "mathusla40": ((0.196, 0.804, 0.196), "dashed (full simulation)",
                   ["Drawn from 10 to 60 GeV; the caption notes the lower masses are under study."]),
    "ship": ((0.855, 0.439, 0.839), "dashed (full simulation)",
             ["The contour is closed by the authors along m_S = 0.021 GeV and sin(theta) = 0.0316,",
              "outside and inside the published frame respectively; kept as published, the",
              "renderer's x range crops the off-frame part. Compared with SHiP's BC4 curve",
              "(ship_bc4.dat) the lower branch is 11x/61x/10x/1.8x deeper at 0.5/1/2/3 GeV:",
              "the quartic also opens B -> K S S. This is expected, not a digitisation error."]),
}


def main() -> int:
    page = fitz.open(str(require(BC5_PDF)))[0]
    for stem, (color, style, notes) in CURVES.items():
        name = {"anubis": "ANUBIS", "codexb": "CODEX-b", "mathusla40": "MATHUSLA40", "ship": "SHiP"}[stem]
        chains = pdfio.merge_chains(pdfio.collect_strokes(page, color, t_samples=SAMPLES, width=2.0), tol=1.5)
        write(PROJ_DIR / f"{stem}_pbc_bc5.dat", chains, pbc_to_data,
              [f"{name} BC5 dark-scalar sensitivity, BR(h -> SS) = 0.01",
               f"Extracted from arXiv:2505.00947v2 / {BC5_SOURCE}",
               f"Date: {DATE}",
               f"Line style in the source: {style}.",
               *notes,
               "Source y axis is sin(theta); values below are squared to sin^2(theta)."])
    return 0


if __name__ == "__main__":
    sys.exit(main())
