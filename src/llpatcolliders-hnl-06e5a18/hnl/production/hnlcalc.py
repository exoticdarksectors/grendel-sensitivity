"""HNLCalc construction helpers."""

from __future__ import annotations

import sys
from pathlib import Path

HNL_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HNL_ROOT / "vendored" / "HNLCalc"))


def init_hnlcalc(flavor):
    from HNLCalc import HNLCalc

    if flavor == "Ue":
        return HNLCalc(ve=1, vmu=0, vtau=0)
    if flavor == "Umu":
        return HNLCalc(ve=0, vmu=1, vtau=0)
    if flavor == "Utau":
        return HNLCalc(ve=0, vmu=0, vtau=1)
    raise ValueError(f"Unknown flavor: {flavor}")
