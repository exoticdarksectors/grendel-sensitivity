"""The HNL scan: how the driver in ``grendel.scan`` is configured for it.

Per (flavor, mass) point the four-vectors are ``<vectors>/<flavor>/combined/
mN_<label>.csv`` and the rest-frame decay templates
``<templates>/<flavor>/templates_<label>.npz`` (exHad, or the FairShip set
behind the first published curve); the template bundle also carries
ctau(|U|^2 = 1), so no separate lifetime table is used.

Policies that fix the published numbers:

* one generator per point, seeded from md5 of ``<flavor>/<label>[/<salt>]``;
* an optional weighted resample of the hitting events (approximate mode),
  drawn from that generator before the acceptance Monte Carlo;
* the event loop may be chunked to bound memory; every chunk draws from the
  same generator, and the yield curve is accumulated chunk by chunk;
* the exclusion band is read off the 200-point log grid
  (``find_exclusion_band``).

The decay-model band co-varies one width nuisance delta(m) through the
lifetime and the hadronic/leptonic template mix on the same built Monte
Carlo (``width_delta``, ``had_frac``); see ``grendel.band.hnl``.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

from ...constants import DEFAULT_FLAVORS, FLAVORS
from ...io.paths import ModelPaths
from ...reco.acceptance import build_event_mc, classify_template_modes
from ...reco.exclusion import find_exclusion_band
from ...scan import (CouplingGrid, ExtraScan, MassPoint, ModelSpec, ScanArrays,
                     ScanConfig, secondary_threshold_columns)

BAND_FIELDS = ("u2_min", "u2_max", "peak_u2", "has_sensitivity",
               "u2_min_open", "u2_max_open")


def seed_for(flavor: str, mass_label: str, seed_salt: str = "") -> int:
    """Deterministic per-point seed (order-independent, reproducible)."""
    key = f"{flavor}/{mass_label}"
    if seed_salt:
        key += f"/{seed_salt}"
    h = hashlib.md5(key.encode()).digest()
    return int.from_bytes(h[:4], "little")


def select_hit_sample(idx, weights, max_hit_events, rng):
    """Hit indices and scan weights for exact or approximate evaluation."""
    weights = np.asarray(weights, dtype=float)
    if max_hit_events is None or max_hit_events <= 0 or len(idx) <= max_hit_events:
        return idx, weights, "exact"

    total_weight = float(weights.sum())
    if total_weight > 0.0 and np.isfinite(total_weight):
        probabilities = weights / total_weight
        sampled_idx = rng.choice(idx, size=int(max_hit_events), replace=True,
                                 p=probabilities)
        sampled_weights = np.full(len(sampled_idx), total_weight / len(sampled_idx),
                                  dtype=float)
        return sampled_idx, sampled_weights, "weighted_resample"

    sampled_pos = rng.choice(len(idx), size=int(max_hit_events), replace=False)
    sampled_idx = idx[sampled_pos]
    sampled_weights = np.asarray(weights[sampled_pos], dtype=float)
    sampled_weights *= len(idx) / max(len(sampled_idx), 1)
    return sampled_idx, sampled_weights, "uniform_rescale"


class TemplateBackend:
    """Decays drawn from a cached rest-frame template bundle."""

    def __init__(self, templates: dict, name: str = "templates"):
        self.templates = templates
        self.name = name
        self.ctau_ref = float(templates["ctau_m_u2eq1"])

    def build(self, p4, direction, entry_d, exit_d, n_samples, rng, return_mc=False):
        out = build_event_mc(p4, direction, entry_d, exit_d, self.templates,
                             n_samples, rng, return_mc=return_mc)
        return out if return_mc else (*out, None)

    def sample_weights(self, template_index):
        return None

    def mode_is_hadronic(self):
        return classify_template_modes(self.templates)


def load_template_bundle(path: Path) -> dict | None:
    """The ``.npz`` bundle as a dict of arrays, or None when absent."""
    path = Path(path)
    if not path.exists():
        return None
    with np.load(path) as archive:
        return {name: archive[name] for name in archive.files}


class HNLSpec(ModelSpec):
    name = "hnl"
    grid = CouplingGrid(-12.0, -1.0, 200)
    scan_mode = "accumulate"
    decay_samples_default = 100
    flavors = tuple(FLAVORS)
    default_flavors = tuple(DEFAULT_FLAVORS)

    def __init__(self, paths: ModelPaths, *, width_delta=None, had_frac=None):
        self.paths = paths
        self.width_delta = width_delta
        self.had_frac = had_frac

    # --- inputs -----------------------------------------------------------
    def vectors_path(self, pt: MassPoint) -> Path:
        return self.paths.vectors / pt.flavor / "combined" / f"mN_{pt.label}.csv"

    def geometry_cache_path(self, pt: MassPoint, vectors_path: Path) -> Path:
        return self.paths.geometry / pt.flavor / f"geom_{pt.label}.npz"

    def templates_path(self, pt: MassPoint) -> Path:
        return self.paths.templates / pt.flavor / f"templates_{pt.label}.npz"

    def decay_backend(self, pt: MassPoint, vectors_path: Path):
        templates = load_template_bundle(self.templates_path(pt))
        if templates is None:
            print(f"  NOTE: no decay templates for {pt.tag}; "
                  "generate them first. Skipping.")
            return None
        backend = TemplateBackend(templates)
        if backend.ctau_ref <= 0:
            return None
        return backend

    # --- policies -----------------------------------------------------------
    def rng_for(self, pt: MassPoint, cfg: ScanConfig):
        return np.random.default_rng(seed_for(pt.flavor, pt.label, cfg.seed_salt))

    def select_events(self, pt, idx, weights, cfg, rng):
        return select_hit_sample(idx, weights, cfg.max_hit_events, rng)

    def chunks(self, pt, n_events, cfg, rng):
        chunk_size = n_events
        if cfg.event_chunk is not None and cfg.event_chunk > 0:
            chunk_size = min(int(cfg.event_chunk), n_events)
        return [(slice(start, min(start + chunk_size, n_events)), rng)
                for start in range(0, n_events, chunk_size)]

    def extra_scans(self, pt, backend, cfg) -> list[ExtraScan]:
        if not self.width_delta:
            return []
        mode_is_had = backend.mode_is_hadronic()
        hf = self.had_frac if self.had_frac else 1.0
        scans = []
        for tag, gscale in (("lo", 1.0 + self.width_delta), ("hi", 1.0 - self.width_delta)):
            w_lep = 1.0 / gscale
            w_had = (1.0 + (gscale - 1.0) / hf) / gscale
            scans.append(ExtraScan(
                tag, gscale,
                lambda tmpl_idx, w_had=w_had, w_lep=w_lep: np.where(
                    mode_is_had[tmpl_idx], w_had, w_lep)))
        return scans

    # --- rows ---------------------------------------------------------------
    def base_row(self, pt, n_events, n_hits):
        return {"mass_GeV": pt.mass, "flavor": pt.flavor,
                "n_events": n_events, "n_hits": n_hits}

    def empty_row(self, pt, n_events, n_hits, cfg, *, evaluated=False):
        row = {"mass_GeV": pt.mass, "flavor": pt.flavor,
               "u2_min": np.nan, "u2_max": np.nan,
               "u2_min_open": False, "u2_max_open": False,
               "peak_N": 0.0, "peak_u2": np.nan,
               "has_sensitivity": False, "n_events": n_events, "n_hits": n_hits}
        if evaluated:
            row.update(n_hits_eval=0, decay_samples=cfg.decay_samples,
                       hit_estimator="exact", event_chunk=cfg.event_chunk,
                       seed_salt=cfg.seed_salt)
        return row

    def finish(self, pt, base, backend, arrays: ScanArrays, cfg: ScanConfig):
        result = find_exclusion_band(arrays.grid, arrays.N, cfg.primary)
        secondary_threshold_columns(
            result, lambda t: find_exclusion_band(arrays.grid, arrays.N, t),
            cfg.secondary, BAND_FIELDS)
        if self.width_delta:
            hf = self.had_frac if self.had_frac else 1.0
            for tag in ("lo", "hi"):
                r_v = find_exclusion_band(arrays.grid, arrays.extras[tag], cfg.primary)
                result[f"u2_max_dm_{tag}"] = r_v["u2_max"]
                result[f"u2_min_dm_{tag}"] = r_v["u2_min"]
            result["width_delta"] = self.width_delta
            result["had_frac"] = hf
        result["mass_GeV"] = pt.mass
        result["flavor"] = pt.flavor
        result["n_events"] = base["n_events"]
        result["n_hits"] = base["n_hits"]
        result["n_hits_eval"] = arrays.n_hits_eval
        result["decay_samples"] = cfg.decay_samples
        result["hit_estimator"] = arrays.hit_estimator
        result["event_chunk"] = cfg.event_chunk
        result["seed_salt"] = cfg.seed_salt
        return result

    # --- optional -----------------------------------------------------------
    def plot(self, csv_path, out_dir):
        from .plot_exclusion import plot_exclusion
        plot_exclusion(csv_path, out_dir)

    def skip_reason(self, pt):
        csv = self.vectors_path(pt)
        if not csv.exists() or csv.stat().st_size == 0:
            return "missing_or_empty_combined_csv"
        if not self.templates_path(pt).exists():
            return "missing_decay_templates"
        return "no_acceptance_or_nonpositive_ctau"
