"""The BC4 scan: how the driver in ``grendel.scan`` is configured for it.

Per mass point the four-vectors are ``<vectors>/mS_<label>.csv`` (inclusive
b -> X_s S over B+, B0, Bs and Lambda_b, ``production.py``). The decay is
either the analytic two-body engine (``acceptance.py``, the first published
curve) or a cached rest-frame template bundle
``<templates>/templates_<label>.npz`` (exHad), whose ctau(sin^2 theta = 1)
then replaces the analytic lifetime so mode composition and lifetime come
from one width model.

Policies that fix the published numbers:

* one generator per mass, seeded ``1000 + i`` for the i-th grid mass and
  ``100000 + round(1000 m)`` off the grid (stable across subset scans); the
  uncertainty campaign passes its own generator;
* no event chunking (the analytic engine bounds memory internally with a
  25 000-decay reconstruction chunk, which is part of its RNG order);
* the exclusion band is refined off the 200-point log grid over
  sin^2 theta in [1e-12, 1e-2] (``find_exclusion_band_refined``).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ...io.paths import ModelPaths
from ...reco.exclusion import find_exclusion_band_refined
from ...scan import (CouplingGrid, MassPoint, ModelSpec, ScanArrays, ScanConfig,
                     secondary_threshold_columns)
from ..hnl.spec import TemplateBackend, load_template_bundle
from . import production
from .acceptance import AnalyticBackend

BAND_FIELDS = ("u2_min", "u2_max", "peak_u2", "has_sensitivity",
               "u2_min_open", "u2_max_open")
DIAGNOSTIC_POINTS = ("u2_min", "peak_u2", "u2_max")


def reconstruction_seed(mass) -> int:
    """Keep canonical mass seeds stable across full-grid and subset scans."""
    rounded = round(float(mass), 3)
    if rounded in production.MASS_GRID:
        return 1000 + production.MASS_GRID.index(rounded)
    return 100_000 + int(round(rounded * 1000))


class ScalarSpec(ModelSpec):
    name = "bc4"
    grid = CouplingGrid(-12.0, -2.0, 200)
    scan_mode = "concatenate"
    decay_samples_default = 100
    flavors = None

    def __init__(self, paths: ModelPaths, *, templates_dir=None, width_scheme="winkler"):
        self.paths = paths
        self.templates_dir = Path(templates_dir) if templates_dir else None
        self.width_scheme = width_scheme

    # --- inputs -----------------------------------------------------------
    def vectors_path(self, pt: MassPoint) -> Path:
        return self.paths.vectors / f"mS_{pt.label}.csv"

    def geometry_cache_path(self, pt: MassPoint, vectors_path: Path) -> Path:
        return self.paths.geometry / f"geom_mS_{pt.label}.npz"

    def decay_backend(self, pt: MassPoint, vectors_path: Path):
        if self.templates_dir is None:
            return AnalyticBackend(pt.mass, self.width_scheme)
        templates = load_template_bundle(self.templates_dir / f"templates_{pt.label}.npz")
        if templates is None:
            print(f"  NOTE: no decay templates for m_S={pt.mass:.3f} in "
                  f"{self.templates_dir}; skipping.", flush=True)
            return None
        backend = TemplateBackend(templates, str(templates.get("decay_model", "templates")))
        if not np.isfinite(backend.ctau_ref) or backend.ctau_ref <= 0:
            return None
        return backend

    # --- policies -----------------------------------------------------------
    def rng_for(self, pt: MassPoint, cfg: ScanConfig):
        return np.random.default_rng(reconstruction_seed(pt.mass) + cfg.seed_offset)

    # --- rows ---------------------------------------------------------------
    def base_row(self, pt, n_events, n_hits):
        return {"mass_GeV": pt.mass, "n_events": int(n_events), "n_hits": n_hits}

    def empty_row(self, pt, n_events, n_hits, cfg, *, evaluated=False):
        return {**self.base_row(pt, n_events, n_hits),
                "u2_min": np.nan, "u2_max": np.nan,
                "u2_min_open": False, "u2_max_open": False,
                "peak_N": 0.0, "peak_u2": np.nan, "has_sensitivity": False}

    def finish(self, pt, base, backend, arrays: ScanArrays, cfg: ScanConfig):
        result = find_exclusion_band_refined(arrays.grid, arrays.N, arrays.evaluate, cfg.primary)
        secondary_threshold_columns(
            result,
            lambda t: find_exclusion_band_refined(arrays.grid, arrays.N, arrays.evaluate, t),
            cfg.secondary, BAND_FIELDS)
        result.update(base)
        result["ctau_sin2th1_m"] = backend.ctau_ref
        result["decay_backend"] = backend.name
        for label in DIAGNOSTIC_POINTS:
            coupling = result[label]
            if not np.isfinite(coupling):
                for field in ("sample_ess", "event_ess", "max_event_fraction"):
                    result[f"{label}_{field}"] = np.nan
                continue
            diagnostics = arrays.diagnostics(coupling)
            for field in ("sample_ess", "event_ess", "max_event_fraction"):
                result[f"{label}_{field}"] = diagnostics[field]
        return result

    # --- optional -----------------------------------------------------------
    def plot(self, csv_path, out_dir):
        from .plot_exclusion import plot_island
        plot_island(csv_path, out_dir)

    def skip_reason(self, pt):
        csv = self.vectors_path(pt)
        if not csv.exists() or csv.stat().st_size == 0:
            return "missing_or_empty_vectors"
        if self.templates_dir is not None and \
                not (self.templates_dir / f"templates_{pt.label}.npz").exists():
            return "missing_decay_templates"
        return "no_acceptance_or_nonpositive_ctau"
