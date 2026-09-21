"""BC10 (fermiophilic ALP) reference curves.

Conventions: ``1/f`` is stored in the BNT normalisation of the ALP model
(g_aff = c_f m_f / f, arXiv:1708.00443)::

    1/f = y / v_h      with  y = 2 v_h / f_GKOZ  and  v_h = 246 GeV
    1/f = g_Y / v_h    (ALPINIST g_Y == y; the ALPINIST LHCb csv over the
                        GKOZ extraction below gives 232 +- 23 across
                        0.35-2.9 GeV, i.e. v_h within digitisation noise)

Existing exclusion (LHCb)
-------------------------
arXiv:2505.00947 has no BC10 figure and GKOZ (arXiv:2310.03524, the paper
behind the model's widths) provide one existing bound: their recast of the
model-independent LHCb searches B -> K + (FIP -> mu mu) for the
fermion-coupled ALP, ``LHCb-constraints.pdf``. The dark-blue solid contour
is the current phenomenology (red = outdated leptons-only, light lines =
decay-length contours). The boundary is clipped at the figure's top frame
(2 v_h/f = 0.8); every chain starts/ends there except the first, which
enters through the left edge (m_a = 0.3 GeV). The excluded region lies
above each chain; the gaps along the top (~3.0-3.2, ~3.6-3.9 GeV) are the
LHCb charmonium vetoes.

Projections and past beam dumps
-------------------------------
Vector-extracted from Ovchynnikov-Zaporozhchenko arXiv:2501.04525 (PRD 112,
015001), ``parameter-space-ALP-fermion.pdf`` (m_a 0.6-2.0 GeV, y^2 axis):
SHiP (blue solid; the dashed old phenomenology is skipped), DarkQuest
phase I (green solid), past beam dumps (black solid, an EXISTING
exclusion) and the LHCb boundary with the revised phenomenology (gray, a
cross-check on the GKOZ extraction). No NA62-dump projection exists with
the post-GKOZ phenomenology (ALPINIST ships dark-scalar contours only,
SensCalc requires Mathematica), so none is drawn.

The x calibration here comes from tick-label text (0.130 pt = +0.40 % in
1/f); three orders below plot scale, but the same defect ``vector_bc4``
corrects. Fix it if the curve is ever regenerated.

Existing exclusion (NA62, kaon regime)
--------------------------------------
The ALPINIST NA62 K+ -> pi+ + inv bound band (``NA62_{upper,lower}limit_gY.csv``
from ``Figures/Bound_data/gY``, commit in ``sources/alpinist/COMMIT_SHA.txt``).

    python -m curves.tools.vector_bc10
"""
from __future__ import annotations

import sys

import fitz
import numpy as np

from ..plot.data import BC10_DATA
from . import pdfio
from .sources import ALPINIST_DIR, GKOZ_LHCB_PDF, OZ_ALP_PDF, SOURCES, require

PROJ_DIR = BC10_DATA / "projections"
VEC_DIR = BC10_DATA / "vector"
V_H = 246.0
LHCB_DATE = "2026-07-10"
DATE = "2026-07-11"
SAMPLES = pdfio.linspace_samples(7)


# ------------------------------------------------------------------ GKOZ --

# Calibration from the figure's own ticks: x tick 84.945 = 0.5 GeV with
# 88.556 pt/GeV; y tick 43.727 = 2 v_h/f = 1e-1 with 48.2198 pt/decade.
X_HALF_GEV, X_PER_GEV = 84.945, 88.556
Y_DEC_M1, Y_PER_DEC = 43.727, 48.2198


def gkoz_to_data(x, y):
    m = 0.5 + (x - X_HALF_GEV) / X_PER_GEV
    log_2vhf = -1.0 - (y - Y_DEC_M1) / Y_PER_DEC
    return m, 10 ** log_2vhf / V_H


def lhcb_excluded():
    page = fitz.open(str(require(GKOZ_LHCB_PDF)))[0]
    chains = []
    for d in page.get_drawings():
        if d.get("type") != "s" or not d.get("color") or not pdfio.near(d["color"], (0.0, 0.0, 1.0)):
            continue
        if len(d.get("items", [])) <= 5:  # legend sample line
            continue
        pts = pdfio.flatten(d["items"])
        if pts:
            chains.append(pts)
    chains.sort(key=lambda c: min(p[0] for p in c))
    print(f"{len(chains)} chains, sizes: {[len(c) for c in chains]}")
    pdfio.write_chains(
        VEC_DIR / "lhcb_bkmumu_excluded.dat", [[gkoz_to_data(x, y) for x, y in ch] for ch in chains],
        ["LHCb B -> K (a ->) mu mu excluded-region boundary, GKOZ re-interpretation",
         "Extracted from arXiv:2310.03524v3 / LHCb-constraints.pdf (dark-blue solid contour)",
         f"Date: {LHCB_DATE}",
         "Source y axis is 2 v_h/f_GKOZ (v_h = 246 GeV); converted here to the",
         "BC10 package convention 1/f_BNT = 2/f_GKOZ, i.e. invf = y_source/246.",
         "The excluded region lies above each chain; chains are clipped at the",
         "source frame top (invf = 3.25e-3) and left edge (m_a = 0.3 GeV).",
         "Gaps along the top edge are the LHCb charmonium veto windows.",
         'Format: chains separated by blank line; each row is "mass_GeV  invf_GeV_inv"'],
        precision=5, numbered=True, columns="mass_GeV  invf_GeV_inv")


# ------------------------------------------------- Ovchynnikov-Zaporozhchenko --

# x LINEAR: '0.6' label centre 81.25, '2.0' centre 503.0
CX06, CX20 = 81.25, 503.0
PT_PER_GEV = (CX20 - CX06) / 1.4
# y log: exponent rows: 10^-7 at 52.65 ... 10^-13 at 276.15
CY_M7, CY_M13 = 52.65, 276.15
PT_PER_DEC = (CY_M13 - CY_M7) / 6.0


def oz_to_data(x, y):
    m = 0.6 + (x - CX06) / PT_PER_GEV
    y2 = 10 ** (-7.0 - (y - CY_M7) / PT_PER_DEC)
    return m, np.sqrt(y2) / V_H          # (m_a, 1/f BNT)


def drop_frame(chains):
    """Axis frame, tick and legend-sample chains are 2-point segments; the
    physics curves have hundreds of vertices."""
    return [ch for ch in chains if len(ch) > 25]


OZ_HEADER = [
    "Extracted from arXiv:2501.04525 (Ovchynnikov-Zaporozhchenko, PRD 112 015001),",
    "parameter-space-ALP-fermion.pdf (solid = revised phenomenology).",
    f"Date: {DATE}",
    "Source y axis is y^2 = (2 v_h/f_GKOZ)^2; converted to 1/f (BNT) = y/v_h,",
    "v_h = 246 GeV. The visible figure frame covers m_a in [0.6, 2.0] GeV;",
    "path data beneath the PDF clip (m_a < 0.6 or > 2.0) is kept as drawn but",
    "was not visually published -- treat it accordingly.",
]


def oz_curves():
    page = fitz.open(str(require(OZ_ALP_PDF)))[0]

    def chains(color):
        # 10 %-opacity strokes are the uncertainty-band edges, not the curves.
        return drop_frame(pdfio.merge_chains(
            pdfio.collect_strokes(page, color, t_samples=SAMPLES, dashed=False, opaque_only=True), tol=1.0))

    def write(path, color, title):
        pdfio.write_chains(path, [[oz_to_data(x, y) for x, y in ch] for ch in chains(color)], title + OZ_HEADER,
                           columns='Format: chains separated by blank line; each row "mass_GeV  invf_GeV^-1"')

    write(PROJ_DIR / "ship_alp2.dat", (0.0, 0.0, 1.0),
          ["SHiP fermiophilic-ALP (BC10) sensitivity, revised phenomenology"])
    write(PROJ_DIR / "darkquest_alp2.dat", (0.0, 0.667, 0.0),
          ["DarkQuest phase I fermiophilic-ALP (BC10) sensitivity, revised phenomenology"])
    write(VEC_DIR / "beamdumps_past_alp2.dat", (0.0, 0.0, 0.0),
          ["Past proton beam-dump exclusions (CHARM et al.) for the fermiophilic ALP,",
           "revised phenomenology -- EXISTING bound"])
    write(VEC_DIR / "lhcb_alp2_crosscheck.dat", (0.5, 0.5, 0.5),
          ["LHCb B -> K mu mu exclusion boundary with the REVISED phenomenology --",
           "cross-check against vector/lhcb_bkmumu_excluded.dat (GKOZ 2023, which",
           "matches the phenomenology of the GRENDEL ALP model)"])


# ---------------------------------------------------------------- ALPINIST --

def na62_band():
    alpinist = SOURCES / ALPINIST_DIR
    up = np.loadtxt(str(alpinist / "NA62_upperlimit_gY.csv"), delimiter=",")
    lo = np.loadtxt(str(alpinist / "NA62_lowerlimit_gY.csv"), delimiter=",")
    VEC_DIR.mkdir(parents=True, exist_ok=True)
    path = VEC_DIR / "na62_kpia_excluded.dat"
    with path.open("w") as fh:
        fh.write("# NA62 K+ -> pi+ + invisible bound on the fermiophilic ALP (existing)\n")
        fh.write("# From ALPINIST Figures/Bound_data/gY (commit sources/alpinist/COMMIT_SHA.txt),\n")
        fh.write("# axis g_Y (dimensionless) converted to 1/f (BNT) = g_Y / v_h, v_h = 246 GeV.\n")
        fh.write(f"# Date: {DATE}\n")
        fh.write("# Excluded band lies between the two chains (lower, then upper).\n")
        fh.write('# Format: chains separated by blank line; each row "mass_GeV  invf_GeV^-1"\n')
        fh.write(f"# chain: {len(lo)} vertices (lower limit)\n")
        for m, gy in lo:
            fh.write(f"{m:.6e}  {gy / V_H:.6e}\n")
        fh.write("\n")
        fh.write(f"# chain: {len(up)} vertices (upper limit)\n")
        for m, gy in up:
            fh.write(f"{m:.6e}  {gy / V_H:.6e}\n")
    print("wrote", path)


def main() -> int:
    lhcb_excluded()
    oz_curves()
    na62_band()
    return 0


if __name__ == "__main__":
    sys.exit(main())
