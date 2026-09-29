"""The BC5 scan: BC4's policies with the quartic channels and a deeper grid."""
from __future__ import annotations

import numpy as np

from ...constants import L_INT_PB
from ...reco.acceptance import scan_u2
from ...scan import CouplingGrid, ScanArrays, ScanConfig
from . import model
from .spec import ScalarSpec

EDGE_LABELS = ("u2_min", "peak_u2", "u2_max")


def channel_yields(arrays: ScanArrays, u2: float) -> dict:
    """N_signal at ``u2`` from the mixing rows and from the quartic rows."""
    power = arrays.coupling_power
    if power is None:
        power = np.ones(len(arrays.weights))
    out = {}
    for label, mask in (("mixing", power == 1.0), ("quartic", power == 0.0)):
        if not mask.any():
            out[label] = 0.0
            continue
        _, n = scan_u2(arrays.d, arrays.passed, arrays.path, arrays.weights * mask,
                       arrays.beta_gamma, arrays.ctau_ref, L_INT_PB, np.asarray([u2]),
                       sample_w=arrays.sample_w, coupling_power=power)
        out[label] = float(n[0])
    return out


class QuarticScalarSpec(ScalarSpec):
    name = "bc5"
    grid = CouplingGrid(-20.0, -2.0, 360)

    def __init__(self, paths, *, templates_dir=None, width_scheme="winkler",
                 br_hss=model.BR_HSS_BC5):
        super().__init__(paths, templates_dir=templates_dir, width_scheme=width_scheme)
        self.br_hss = br_hss

    def base_row(self, pt, n_events, n_hits):
        return {**super().base_row(pt, n_events, n_hits),
                "br_hss": self.br_hss, "alpha_GeV": model.alpha_quartic(pt.mass, self.br_hss)}

    def finish(self, pt, base, backend, arrays: ScanArrays, cfg: ScanConfig):
        result = super().finish(pt, base, backend, arrays, cfg)
        power = arrays.coupling_power
        result["n_hits_quartic"] = int((power == 0.0).sum()) if power is not None else 0
        for label in EDGE_LABELS:
            u2 = result.get(label, np.nan)
            split = channel_yields(arrays, u2) if np.isfinite(u2) else {"mixing": np.nan, "quartic": np.nan}
            result[f"{label}_N_mixing"] = split["mixing"]
            result[f"{label}_N_quartic"] = split["quartic"]
        return result
