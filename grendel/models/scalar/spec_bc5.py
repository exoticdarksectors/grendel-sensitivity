"""The BC5 scan: BC4's policies with the quartic channels and a deeper grid.

Besides the nominal island the scan solves, on the same Monte Carlo, the island of every production
mode alone and the islands rescaled to other values of BR(h -> SS). Every such curve is a fixed linear
combination of the per-mode yield curves: the mixing rows do not depend on BR(h -> SS), and every
quartic row's weight is linear in it (sigma_h BR for h -> S S, alpha^2 for the quartic B decays).
"""
from __future__ import annotations

import numpy as np

from ...constants import L_INT_PB
from ...reco.acceptance import scan_u2
from ...reco.exclusion import find_exclusion_band_refined
from ...scan import CouplingGrid, ScanArrays, ScanConfig, secondary_threshold_columns
from . import model
from .production_bc5 import CHANNEL_ID, CHANNELS
from .spec import BAND_FIELDS, ScalarSpec

EDGE_LABELS = ("u2_min", "peak_u2", "u2_max")
MODES = CHANNELS
ISLAND_FIELDS = ("u2_min", "u2_max", "peak_N", "peak_u2", "has_sensitivity",
                 "u2_min_open", "u2_max_open")
BR_HSS_OVERLAYS_DEFAULT = (0.001,)


def br_tag(br_hss: float) -> str:
    """Column prefix of the curves rescaled to BR(h -> SS) = ``br_hss``: ``brhss0.001``."""
    return f"brhss{br_hss:g}"


def mode_masks(arrays: ScanArrays) -> dict[str, np.ndarray]:
    """Which hitting rows each production mode owns. The channel column tells ``mixing``, ``hSS`` and
    ``BSS`` apart; without it (CSVs written before the column existed) the coupling power only tells
    ``mixing`` from ``quartic`` = hSS + BSS."""
    if arrays.channel is not None:
        return {name: arrays.channel == CHANNEL_ID[name] for name in MODES}
    power = arrays.coupling_power
    if power is None:
        power = np.ones(len(arrays.weights))
    return {"mixing": power == 1.0, "quartic": power == 0.0}


class ModeYields:
    """N_signal(sin^2 theta) of each production mode, on the scan grid and at any coupling; all the
    modes come out of one pass over the Monte Carlo."""

    def __init__(self, arrays: ScanArrays):
        self.arrays = arrays
        self.masks = mode_masks(arrays)
        self.names = list(self.masks)
        self.weights = np.stack([arrays.weights * self.masks[name] for name in self.names])
        self.on_grid = dict(zip(self.names, self._scan(arrays.grid)))

    def _scan(self, u2):
        """``(mode, u2)`` yields."""
        a = self.arrays
        _, n = scan_u2(a.d, a.passed, a.path, self.weights, a.beta_gamma, a.ctau_ref, L_INT_PB,
                       np.atleast_1d(np.asarray(u2, float)), sample_w=a.sample_w,
                       coupling_power=a.coupling_power)
        return n

    def at(self, u2: float) -> dict[str, float]:
        """N_signal at ``u2`` per mode."""
        return dict(zip(self.names, self._scan(u2)[:, 0].tolist()))

    def combine(self, coefficients: dict[str, float]):
        """``(N on the grid, evaluate)`` of ``sum_mode c_mode N_mode``, as the island solver wants."""
        N = sum((c * self.on_grid[name] for name, c in coefficients.items()),
                np.zeros(len(self.arrays.grid)))

        def evaluate(u2):
            yields = self.at(u2)
            return sum(c * yields[name] for name, c in coefficients.items())

        return N, evaluate


def variant_coefficients(modes, br_hss: float, overlays) -> dict[str, dict[str, float]]:
    """Column prefix -> {mode: coefficient} of every curve solved besides the nominal one: each
    production mode alone at the nominal BR(h -> SS), and at every overlay BR the full curve and each
    quartic mode alone (the mixing mode does not depend on BR, so it is not repeated)."""
    modes = list(modes)
    out = {name: {name: 1.0} for name in modes}
    for br in overlays:
        ratio = br / br_hss
        tag = br_tag(br)
        out[tag] = {name: (1.0 if name == "mixing" else ratio) for name in modes}
        for name in modes:
            if name != "mixing":
                out[f"{tag}_{name}"] = {name: ratio}
    return out


def edge_yields(yields: ModeYields, u2: float) -> dict[str, float]:
    """N_signal at ``u2`` from the mixing rows, the quartic rows and (channel column present) each
    quartic mode; NaN where unknown."""
    out = {name: np.nan for name in ("mixing", "quartic", "hSS", "BSS")}
    if not np.isfinite(u2):
        return out
    out.update(yields.at(u2))
    if "quartic" not in yields.masks:
        out["quartic"] = out["hSS"] + out["BSS"]
    return out


def channel_yields(arrays: ScanArrays, u2: float) -> dict:
    """N_signal at ``u2`` from the mixing rows, the quartic rows and each quartic mode."""
    return edge_yields(ModeYields(arrays), u2)


class QuarticScalarSpec(ScalarSpec):
    name = "bc5"
    grid = CouplingGrid(-20.0, -2.0, 360)

    def __init__(self, paths, *, templates_dir=None, width_scheme="winkler",
                 br_hss=model.BR_HSS_BC5, br_hss_overlays=BR_HSS_OVERLAYS_DEFAULT):
        super().__init__(paths, templates_dir=templates_dir, width_scheme=width_scheme)
        self.br_hss = br_hss
        overlays = tuple(float(br) for br in br_hss_overlays)
        if any(br <= 0.0 for br in overlays):
            raise ValueError(f"BR(h -> SS) overlays must be positive, got {overlays}")
        self.br_hss_overlays = overlays if br_hss > 0.0 else ()

    def base_row(self, pt, n_events, n_hits):
        return {**super().base_row(pt, n_events, n_hits),
                "br_hss": self.br_hss, "alpha_GeV": model.alpha_quartic(pt.mass, self.br_hss)}

    def finish(self, pt, base, backend, arrays: ScanArrays, cfg: ScanConfig):
        result = super().finish(pt, base, backend, arrays, cfg)
        yields = ModeYields(arrays)
        power = arrays.coupling_power
        result["n_hits_quartic"] = int((power == 0.0).sum()) if power is not None else 0
        for name in MODES:
            result[f"n_hits_{name}"] = (int(yields.masks[name].sum())
                                        if name in yields.masks else np.nan)
        for label in EDGE_LABELS:
            for name, value in edge_yields(yields, result.get(label, np.nan)).items():
                result[f"{label}_N_{name}"] = value
        variants = variant_coefficients(yields.masks, self.br_hss, self.br_hss_overlays)
        for prefix, coefficients in variants.items():
            N, evaluate = yields.combine(coefficients)
            island = find_exclusion_band_refined(arrays.grid, N, evaluate, cfg.primary)
            for field in ISLAND_FIELDS:
                result[f"{prefix}_{field}"] = island[field]
            secondary_threshold_columns(
                result, lambda t: find_exclusion_band_refined(arrays.grid, N, evaluate, t),
                cfg.secondary, BAND_FIELDS, rename=lambda src: f"{prefix}_{src}")
        return result
