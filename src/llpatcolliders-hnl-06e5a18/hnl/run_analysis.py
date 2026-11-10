#!/usr/bin/env python3
"""
Top-level entry point for the HNL analysis pipeline.

Mirrors hnl/run_all.py (the production driver). Reads combined HNL 4-vector
CSVs from ``hnl/tmp/runs/<tag>/llp_4vectors/`` (or the directory pointed at
by ``$HNL_LLP_VECTORS_DIR``), computes the (m_N, U^2) sensitivity, and
writes the money plot under ``hnl/tmp/runs/<tag>/analysis/``.

Default flavor is Umu (HNL pattern 010 -- muon-mixing scenario).

Usage:

    python run_analysis.py                       # Umu, full grid, 3 workers
    python run_analysis.py --workers 4           # custom workers
    python run_analysis.py --mass 0.5 1.0 2.0    # subset of masses
    python run_analysis.py --flavor Ue Umu       # multi-flavor scan
    python run_analysis.py --plot-only           # re-plot only
"""

from __future__ import annotations

import sys
from pathlib import Path

HNL_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(HNL_ROOT))
sys.path.insert(0, str(HNL_ROOT / "vendored" / "HNLCalc"))

from analysis.run_sensitivity import main


if __name__ == "__main__":
    sys.exit(main())
