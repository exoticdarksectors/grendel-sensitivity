"""The BC10 scan: how the driver in ``grendel.scan`` is configured for it.

BC10 is coupling-controlled like the HNL: the single inverse decay constant
1/f fixes the production yield, the lifetime and the visible branching
ratios at once. The scan variable is

    u2 = ((1/f) / (1/f_ref))^2          (1/f_ref = model.INV_F_REF)

so that off the reference point at which the four-vectors are weighted and
the templates store ctau, the production yield scales as u2, ctau as 1/u2,
and the visible branching ratios not at all -- the structure ``scan_u2``
reweights. Island edges are mapped back to 1/f for the result rows
(``invf_min`` from ``u2_min``, ``invf_max`` from ``u2_max``). Masses inside
the eta and eta' poles are not scanned; their rows carry ``exclusion_reason``.

Policies that fix the published numbers:

* the acceptance Monte Carlo runs in chunks of 256 hitting events, each with
  its own generator seeded ``round(1000 m) * 100000 + offset + chunk``;
* events are weighted by the visible fraction unless the template bundle
  already samples the full branching mixture; templates may carry per-sample
  matrix-element weights;
* the exclusion band is refined off a 300-point log grid over
  u2 in [1e-16, 1] (``find_exclusion_band_refined``).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ...io.paths import ModelPaths
from ...reco.acceptance import build_event_mc
from ...reco.exclusion import find_exclusion_band_refined
from ...scan import CouplingGrid, MassPoint, ModelSpec, ScanArrays, ScanConfig
from ..hnl.spec import load_template_bundle
from . import model

EVENT_CHUNK = 256
DIAGNOSTIC_POINTS = (("invf_min", "u2_min"), ("peak", "peak_u2"), ("invf_max", "u2_max"))


def invf_from_u2(u2):
    """Map the scan variable back to the physical coupling 1/f [GeV^-1]."""
    return model.INV_F_REF * np.sqrt(u2)


def _invf_or_nan(u2):
    return invf_from_u2(u2) if np.isfinite(u2) else np.nan


class ALPTemplateBackend:
    """Decays drawn from a cached template bundle, with the bundle's
    optional per-sample matrix-element weights."""

    def __init__(self, templates: dict):
        self.templates = templates
        self.name = str(templates["decay_backend"]) if "decay_backend" in templates else "legacy"
        self.ctau_ref = float(templates["ctau_m_u2eq1"])
        self.includes_full_branching = ("includes_full_branching" in templates
                                        and bool(templates["includes_full_branching"]))
        self._matrix_element_weight = (np.asarray(templates["matrix_element_weight"])
                                       if "matrix_element_weight" in templates else None)

    def build(self, p4, direction, entry_d, exit_d, n_samples, rng, return_mc=False):
        out = build_event_mc(p4, direction, entry_d, exit_d, self.templates,
                             n_samples, rng, return_mc=return_mc)
        return out if return_mc else (*out, None)

    def sample_weights(self, template_index):
        if self._matrix_element_weight is None:
            return None
        return self._matrix_element_weight[template_index]


class ALPSpec(ModelSpec):
    name = "bc10"
    grid = CouplingGrid(-16.0, 0.0, 300)
    scan_mode = "concatenate"
    decay_samples_default = 60
    flavors = None

    def __init__(self, paths: ModelPaths):
        self.paths = paths

    # --- inputs -----------------------------------------------------------
    def pre_check(self, pt: MassPoint):
        resonance = model.excluded_light_meson_resonance(pt.mass)
        if resonance is None:
            return None
        return {"mass_GeV": pt.mass, "n_events": 0, "n_hits": 0,
                "has_sensitivity": False, "peak_N": np.nan,
                "invf_min": np.nan, "invf_max": np.nan,
                "invf_min_open": False, "invf_max_open": False,
                "peak_invf": np.nan,
                "exclusion_reason": f"unsupported {resonance} resonance"}

    def vectors_path(self, pt: MassPoint) -> Path:
        return self.paths.vectors / f"mA_{pt.label}.csv"

    def geometry_cache_path(self, pt: MassPoint, vectors_path: Path) -> Path:
        return self.paths.geometry / f"geom_{pt.label}.npz"

    def templates_path(self, pt: MassPoint) -> Path:
        return self.paths.templates / f"templates_{pt.label}.npz"

    def decay_backend(self, pt: MassPoint, vectors_path: Path):
        templates = load_template_bundle(self.templates_path(pt))
        if templates is None:
            print(f"  m_a={pt.mass:.3f}: missing templates; generate them first")
            return None
        backend = ALPTemplateBackend(templates)
        if backend.ctau_ref <= 0:
            return None
        return backend

    # --- policies -----------------------------------------------------------
    def event_weights(self, pt, backend, weights):
        factor = 1.0 if backend.includes_full_branching else model.visible_fraction(pt.mass)
        return weights * factor

    def chunks(self, pt, n_events, cfg, rng):
        base_seed = int(round(pt.mass * 1000)) * 100_000 + int(cfg.seed_offset)
        return [(slice(cs, min(cs + EVENT_CHUNK, n_events)), np.random.default_rng(base_seed + ci))
                for ci, cs in enumerate(range(0, n_events, EVENT_CHUNK))]

    # --- rows ---------------------------------------------------------------
    def base_row(self, pt, n_events, n_hits):
        return {"mass_GeV": pt.mass, "n_events": n_events, "n_hits": n_hits}

    def empty_row(self, pt, n_events, n_hits, cfg, *, evaluated=False):
        return {**self.base_row(pt, n_events, n_hits), "has_sensitivity": False,
                "peak_N": 0.0, "invf_min": np.nan, "invf_max": np.nan,
                "invf_min_open": False, "invf_max_open": False, "peak_invf": np.nan}

    def finish(self, pt, base, backend, arrays: ScanArrays, cfg: ScanConfig):
        band = find_exclusion_band_refined(arrays.grid, arrays.N, arrays.evaluate, cfg.primary)
        res = {
            **base,
            "decay_backend": backend.name,
            "has_sensitivity": band["has_sensitivity"],
            "peak_N": band["peak_N"],
            "peak_invf": invf_from_u2(band["peak_u2"]),
            "invf_min": _invf_or_nan(band["u2_min"]),
            "invf_max": _invf_or_nan(band["u2_max"]),
            "invf_min_open": bool(band["u2_min_open"]),
            "invf_max_open": bool(band["u2_max_open"]),
        }
        for threshold in cfg.secondary:
            tag = f"N{threshold:g}"
            extra = find_exclusion_band_refined(arrays.grid, arrays.N, arrays.evaluate, threshold)
            res[f"has_sensitivity_{tag}"] = extra["has_sensitivity"]
            res[f"peak_invf_{tag}"] = invf_from_u2(extra["peak_u2"])
            res[f"invf_min_{tag}"] = _invf_or_nan(extra["u2_min"])
            res[f"invf_max_{tag}"] = _invf_or_nan(extra["u2_max"])
            res[f"invf_min_open_{tag}"] = bool(extra["u2_min_open"])
            res[f"invf_max_open_{tag}"] = bool(extra["u2_max_open"])
        for label, key in DIAGNOSTIC_POINTS:
            u2 = band[key]
            if not np.isfinite(u2):
                for field in ("sample_ess", "event_ess", "max_event_fraction"):
                    res[f"{label}_{field}"] = np.nan
                continue
            diagnostics = arrays.diagnostics(u2)
            for field in ("sample_ess", "event_ess", "max_event_fraction"):
                res[f"{label}_{field}"] = diagnostics[field]
        return res

    # --- optional -----------------------------------------------------------
    def plot(self, csv_path, out_dir):
        from .plot import plot_island
        plot_island(csv_path, out_dir)

    def skip_reason(self, pt):
        csv = self.vectors_path(pt)
        if not csv.exists() or csv.stat().st_size == 0:
            return "missing_or_empty_vectors"
        if not self.templates_path(pt).exists():
            return "missing_decay_templates"
        return "no_acceptance_or_nonpositive_ctau"
