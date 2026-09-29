"""The BC4 scan: how the driver in ``grendel.scan`` is configured for it."""
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

    def vectors_path(self, pt: MassPoint) -> Path:
        return self.paths.vectors / f"mS_{pt.label}.csv"

    def geometry_cache_path(self, pt: MassPoint, vectors_path: Path) -> Path:
        return self.paths.geometry / f"geom_mS_{pt.label}.npz"

    def input_files(self, pt: MassPoint) -> list[Path]:
        files = [self.vectors_path(pt)]
        if self.templates_dir is not None:
            files.append(self.templates_dir / f"templates_{pt.label}.npz")
        return files

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

    def rng_for(self, pt: MassPoint, cfg: ScanConfig):
        return np.random.default_rng(reconstruction_seed(pt.mass) + cfg.seed_offset)

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

    def skip_reason(self, pt):
        csv = self.vectors_path(pt)
        if not csv.exists() or csv.stat().st_size == 0:
            return "missing_or_empty_vectors"
        if self.templates_dir is not None and \
                not (self.templates_dir / f"templates_{pt.label}.npz").exists():
            return "missing_decay_templates"
        return "no_acceptance_or_nonpositive_ctau"
