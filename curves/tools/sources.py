"""The published figures the reference curves are extracted from.

They live under ``curves/sources`` (or ``GRENDEL_CURVES_SOURCES``) at the
paths below, which mirror each arXiv source bundle. None is tracked; the
README in that directory says where to download them.
"""
from __future__ import annotations

import os
from pathlib import Path

SOURCES = Path(os.environ.get("GRENDEL_CURVES_SOURCES",
                              Path(__file__).resolve().parents[1] / "sources")).expanduser()

PBC_2025_FIGURES = "PBC/arXiv-2505.00947v2/Figures/5-PhysicsReach/BSM_benchmarks"
BC7_PDF = f"{PBC_2025_FIGURES}/pbc_bc7_gray_final.pdf"
BC5_PDF = f"{PBC_2025_FIGURES}/pbc_bc5_gray_final.pdf"
GKOZ_LHCB_PDF = "arXiv-2310.03524v3/LHCb-constraints.pdf"
SHIP_BC4_PDF = "arXiv-2504.06692/BC4.pdf"
CODEXB_BDECAY_PDF = "arXiv-1911.00481/B_decay_benchmark.pdf"
CODEXB_HNL_PDFS = {"Ue": "arXiv-1911.00481/HNL_Ele.pdf", "Umu": "arXiv-1911.00481/HNL_Mu.pdf",
                   "Utau": "arXiv-1911.00481/HNL_Tau.pdf"}
OZ_ALP_PDF = "arXiv-2501.04525/parameter-space-ALP-fermion.pdf"
FIPS2022_ALP_PDF = "arXiv-2305.01715v1/Light_DM/summary_plots_LDM/ALPs_fermions_fips2022_text_logo.pdf"
ALPINIST_DIR = "alpinist"


def require(relative: str) -> Path:
    path = SOURCES / relative
    if not path.is_file():
        raise FileNotFoundError(f"source figure not found: {path}\n"
                                f"See curves/sources/README.md, or set GRENDEL_CURVES_SOURCES.")
    return path
