"""BC10 existing exclusions of the FIPs 2022 summary figure.

``ALPs_fermions_fips2022_text_logo.pdf`` of arXiv:2305.01715 (g_Y = 2v/f
vs m_ALP, page rotated 90 degrees) is the most recent comprehensive
compilation of what is already excluded for the fermiophilic ALP. It
predates the GKOZ phenomenology revision (arXiv:2310.03524,
arXiv:2501.04525), so regions that depend on the ALP lifetime or hadronic
widths carry that caveat.

Region groups, identified by fill colour::

    gray (0.8,0.8,0.8) opacity 0.4   K+ -> pi+ X, K+ -> pi+ inv (E949), CHARM,
                                     KL -> pi0 ll (several polygons, split by
                                     bounding box into named regions)
    (0.6,0.6,0.6) opacity 0.2        K+ -> pi+ + inv (E949 data), two lobes
    (1.0,0.4,1.0) opacity 0.3        LHCb B+ -> K+ mu mu
    (0.8,0.4,0.8) opacity 0.3        LHCb B0 -> K*0 mu mu
    (0.6,0.6,0.8)                    BaBar

Output: ``curves/data/bc10/vector/fips2022/<group>.dat`` with rows
``mass_GeV  invf_GeV^-1`` (1/f BNT = g_Y / v_h); g_Y = invf * 246.

    python -m curves.tools.vector_bc10_fips2022
"""
from __future__ import annotations

import sys

import fitz
import numpy as np

from ..plot.data import BC10_DATA
from . import pdfio
from .sources import FIPS2022_ALP_PDF, require

OUT_DIR = BC10_DATA / "vector" / "fips2022"
V_H = 246.0
DATE = "2026-07-12"
SAMPLES = pdfio.linspace_samples(6)

DESCRIPTIONS = {
    "kaon_charm_union": "Union of the kaon-sector bounds (K+ -> pi+ X low part, E949 environs) and the CHARM main region",
    "charm_islet": "CHARM high-mass islet (~1 GeV)",
    "kpix_upper": "K+ -> pi+ X upper wedge",
    "b_to_k_inv": "B -> K + invisible",
    "kl_pi0ll": "KL -> pi0 l l band",
    "bs_mumu": "Bs -> mu mu (top right)",
    "bbn": "BBN constraints (lifetime-based cosmology bound; pre-GKOZ lifetimes -- treat the boundary as indicative)",
}


class Frame:
    """The ROOT frame spans exactly m 1e-2..10, g_Y 1e-8..1; its corners in
    display (rotated) coordinates calibrate the axes."""

    def __init__(self, page):
        self.rot = page.rotation_matrix
        frame = None
        for d in page.get_drawings():
            if d["type"] not in ("s", "fs") or not d["color"]:
                continue
            if tuple(round(c, 3) for c in d["color"]) != (0.0, 0.0, 0.0):
                continue
            r = d["rect"]
            if frame is None or r.width * r.height > frame.width * frame.height:
                frame = r
        fx0, fy0 = self.to_display(frame.x0, frame.y0)
        fx1, fy1 = self.to_display(frame.x1, frame.y1)
        self.x0, self.x1 = min(fx0, fx1), max(fx0, fx1)   # m = 1e-2 .. 10 (3 decades)
        self.y0, self.y1 = min(fy0, fy1), max(fy0, fy1)   # gY = 1 .. 1e-8 (8 decades, top = 1)
        print(f"frame display: x [{self.x0:.1f},{self.x1:.1f}] y [{self.y0:.1f},{self.y1:.1f}]")

    def to_display(self, x, y):
        pt = fitz.Point(x, y) * self.rot
        return pt.x, pt.y

    def to_data(self, x, y):
        xd, yd = self.to_display(x, y)
        m = 10 ** (-2.0 + 3.0 * (xd - self.x0) / (self.x1 - self.x0))
        gy = 10 ** (0.0 - 8.0 * (yd - self.y0) / (self.y1 - self.y0))
        return m, gy / V_H


def classify_gray(frame, ch):
    """Map the polygons of the light-gray group to named regions by their
    bounding boxes (pinned once by eye against the source figure)."""
    m = np.array([frame.to_data(x, y)[0] for x, y in ch])
    gy = np.array([frame.to_data(x, y)[1] * V_H for x, y in ch])
    if gy.max() < 5.0e-5:
        return "bbn"
    if m.min() > 3.0:
        return "bs_mumu"
    if 0.20 < m.min() and m.max() < 0.40:
        return "kl_pi0ll"
    if m.max() > 0.95 and m.max() < 1.5:
        return "charm_islet"
    if gy.min() > 2.0e-2:
        return "kpix_upper"
    if gy.min() > 5.0e-3:
        return "b_to_k_inv"
    return "kaon_charm_union"


def write(frame, name, chains, description):
    pdfio.write_chains(
        OUT_DIR / f"{name}.dat", [[frame.to_data(x, y) for x, y in ch] for ch in chains],
        [description,
         "Extracted from arXiv:2305.01715 (FIPs 2022 report), ALPs_fermions_fips2022_text_logo.pdf",
         f"Date: {DATE}",
         "NOTE: pre-GKOZ phenomenology (see module docstring).",
         "Source y axis g_Y = 2 v_h/f; stored as 1/f (BNT) = g_Y/246."],
        columns='Format: chains separated by blank line; each row "mass_GeV  invf_GeV^-1"')


def main() -> int:
    page = fitz.open(str(require(FIPS2022_ALP_PDF)))[0]
    frame = Frame(page)
    groups: dict[str, list] = {}
    for ch in pdfio.collect_fills(page, (0.8, 0.8, 0.8), t_samples=SAMPLES, opacity=0.4, rects=True):
        groups.setdefault(classify_gray(frame, ch), []).append(ch)
    for name, chains in groups.items():
        write(frame, name, chains, DESCRIPTIONS[name] + " -- existing exclusion")
    write(frame, "e949_kpi_inv", pdfio.collect_fills(page, (0.6, 0.6, 0.6), t_samples=SAMPLES, opacity=0.2, rects=True),
          "K+ -> pi+ + inv (based on E949 data), two lobes around the pi+ veto gap -- existing exclusion")
    write(frame, "lhcb_kmumu", pdfio.collect_fills(page, (1.0, 0.4, 1.0), t_samples=SAMPLES, opacity=0.3, rects=True),
          "LHCb B+ -> K+ mu mu excluded (FIPs 2022 version; cf. the GKOZ recast in ../lhcb_bkmumu_excluded.dat)")
    write(frame, "lhcb_kstmumu", pdfio.collect_fills(page, (0.8, 0.4, 0.8), t_samples=SAMPLES, opacity=0.3, rects=True),
          "LHCb B0 -> K*0 mu mu excluded (FIPs 2022 version)")
    write(frame, "babar", pdfio.collect_fills(page, (0.6, 0.6, 0.8), t_samples=SAMPLES, rects=True),
          "BaBar excluded region")
    return 0


if __name__ == "__main__":
    sys.exit(main())
