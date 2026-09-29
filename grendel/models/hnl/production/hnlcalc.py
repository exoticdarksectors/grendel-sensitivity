"""HNLCalc construction helpers."""

from __future__ import annotations

import sys

from ....io.paths import third_party_dir


def _import_hnlcalc():
    path = third_party_dir() / "HNLCalc"
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))
    from HNLCalc import HNLCalc
    return HNLCalc


def init_hnlcalc(flavor):
    HNLCalc = _import_hnlcalc()
    if flavor == "Ue":
        return HNLCalc(ve=1, vmu=0, vtau=0)
    if flavor == "Umu":
        return HNLCalc(ve=0, vmu=1, vtau=0)
    if flavor == "Utau":
        return HNLCalc(ve=0, vmu=0, vtau=1)
    raise ValueError(f"Unknown flavor: {flavor}")
