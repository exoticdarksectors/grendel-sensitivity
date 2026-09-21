"""The per-mass sensitivity scan every benchmark model runs through.

A model supplies a ``ModelSpec``: where its four-vectors and decay templates
are, how it seeds and chunks the acceptance Monte Carlo, which solver it
uses and how it names its columns. The driver owns everything the models
share -- loading, the cached ray-cast, the decay backend, the chunked event
loop and the coupling scan -- and hands the scan arrays back to the model to
turn into a result row.

Seeding, chunking and the choice of solver decide the last digit of every
published number, so they are model policy, deliberately not unified here.

Result layout, per model, under the analysis directory::

    sensitivity.csv            one row per mass point (and flavor)
    sensitivity.partial.csv    checkpoint while a scan runs
    run_metadata.json          configuration, coverage and timing
    scan_status.json           progress
    geometry_cache/            ray-cast cache (unless redirected)
"""
from __future__ import annotations

import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Callable, Iterable, Protocol

import numpy as np
import pandas as pd

from .constants import L_INT_PB
from .geometry.raycast import directions_from_eta_phi, get_mesh, load_or_compute_geometry
from .io.atomic import atomic_csv, atomic_json
from .io.vectors import format_mass_for_filename, load_combined_csv
from .reco.acceptance import P_CUT, scan_u2, signal_contribution_diagnostics


# ---------------------------------------------------------------- inputs --

@dataclass(frozen=True)
class MassPoint:
    mass: float
    flavor: str | None = None

    @property
    def label(self) -> str:
        return format_mass_for_filename(self.mass)

    @property
    def tag(self) -> str:
        return f"{self.flavor}/{self.label}" if self.flavor else self.label

    @property
    def key(self) -> tuple:
        return (self.flavor, float(self.mass))


@dataclass(frozen=True)
class ScanConfig:
    """Knobs a scan takes from the command line."""
    decay_samples: int
    thresholds: tuple[float, ...] = (3.0,)   # primary first; the rest are solved on the same yield curve
    max_hit_events: int | None = None        # weighted resampling cap (HNL approximate mode)
    event_chunk: int | None = None           # memory bound for the event loop, when the model chunks
    seed_salt: str = ""                      # independent numerical-control repeats
    seed_offset: int = 0                     # added to every per-mass seed (BC10 campaigns)
    force_geometry: bool = False
    save_diagnostics: bool = False

    @property
    def primary(self) -> float:
        return self.thresholds[0]

    @property
    def secondary(self) -> tuple[float, ...]:
        return tuple(self.thresholds[1:])


@dataclass(frozen=True)
class CouplingGrid:
    log10_min: float
    log10_max: float
    n: int

    def values(self) -> np.ndarray:
        return np.logspace(self.log10_min, self.log10_max, self.n)


@dataclass(frozen=True)
class ExtraScan:
    """A co-varied yield curve accumulated alongside the nominal one: the
    lifetime scaled by ``ctau_scale`` and each sample reweighted by
    ``sample_weight(template_index)``."""
    tag: str
    ctau_scale: float
    sample_weight: Callable[[np.ndarray], np.ndarray]


class DecayBackend(Protocol):
    """One mass point's decay model: how daughters are drawn, and the
    lifetime at unit coupling the scan reweights from."""
    name: str
    ctau_ref: float

    def build(self, p4, direction, entry_d, exit_d, n_samples, rng, return_mc=False):
        """-> (d, passed, template_index or None, mc or None); ``mc`` holds the
        full-length selection arrays (``acceptance.scatter_mc``) when asked for."""

    def sample_weights(self, template_index):
        """Per-sample weights (matrix-element reweighting) or None."""


# ------------------------------------------------------------ the result --

@dataclass
class ScanArrays:
    """Everything the coupling scan produced for one mass point."""
    grid: np.ndarray
    N: np.ndarray
    d: np.ndarray
    passed: np.ndarray
    path: np.ndarray
    weights: np.ndarray
    beta_gamma: np.ndarray
    ctau_ref: float
    sample_w: np.ndarray | None = None
    extras: dict[str, np.ndarray] = field(default_factory=dict)
    n_hits_eval: int = 0
    hit_estimator: str = "exact"
    n_samples: int = 0
    template_index: np.ndarray | None = None
    mc: dict | None = None          # selection arrays, only with return_mc

    def evaluate(self, u2: float) -> float:
        _, signal = scan_u2(self.d, self.passed, self.path, self.weights,
                            self.beta_gamma, self.ctau_ref, L_INT_PB,
                            np.asarray([u2]), sample_w=self.sample_w)
        return signal[0]

    def diagnostics(self, u2: float) -> dict:
        return signal_contribution_diagnostics(
            self.d, self.passed, self.path, self.weights, self.beta_gamma,
            self.ctau_ref, u2, sample_w=self.sample_w)


@dataclass
class PointResult:
    row: dict
    arrays: ScanArrays | None = None


# ------------------------------------------------------------- the model --

class ModelSpec:
    """What a benchmark model tells the driver. Subclasses override the
    policies; the defaults are the simplest choice (exact hits, one chunk,
    identity weights)."""

    name: str = ""
    grid: CouplingGrid = CouplingGrid(-12.0, -1.0, 200)
    scan_mode: str = "concatenate"   # "concatenate": one scan over all chunks; "accumulate": scan per chunk and sum
    decay_samples_default: int = 100
    flavors: tuple[str, ...] | None = None

    # --- inputs -----------------------------------------------------------
    def points(self, masses, flavors=None) -> list[MassPoint]:
        if self.flavors is None:
            return [MassPoint(float(m)) for m in masses]
        return [MassPoint(float(m), f) for f in (flavors or self.flavors) for m in masses]

    def pre_check(self, pt: MassPoint) -> dict | None:
        """A result row that short-circuits the point, or None."""
        return None

    def vectors_path(self, pt: MassPoint) -> Path:
        raise NotImplementedError

    def geometry_cache_path(self, pt: MassPoint, vectors_path: Path) -> Path:
        raise NotImplementedError

    def decay_backend(self, pt: MassPoint, vectors_path: Path) -> DecayBackend | None:
        """None means the point is skipped (missing templates, closed decay)."""
        raise NotImplementedError

    # --- policies that fix the numbers ---------------------------------------
    def rng_for(self, pt: MassPoint, cfg: ScanConfig):
        """The generator the whole point draws from (None: per-chunk seeding)."""
        return None

    def select_events(self, pt, idx, weights, cfg, rng):
        """-> (idx, weights, estimator): exact hits, or a weighted resample."""
        return idx, np.asarray(weights, dtype=float), "exact"

    def event_weights(self, pt, backend, weights):
        return weights

    def chunks(self, pt, n_events, cfg, rng) -> Iterable[tuple[slice, np.random.Generator]]:
        return [(slice(0, n_events), rng)]

    def extra_scans(self, pt, backend, cfg) -> list[ExtraScan]:
        return []

    # --- rows ---------------------------------------------------------------
    def base_row(self, pt: MassPoint, n_events: int, n_hits: int) -> dict:
        raise NotImplementedError

    def empty_row(self, pt, n_events, n_hits, cfg, *, evaluated=False) -> dict:
        raise NotImplementedError

    def finish(self, pt, base: dict, backend, arrays: ScanArrays, cfg: ScanConfig) -> dict:
        raise NotImplementedError

    # --- optional -----------------------------------------------------------
    def plot(self, csv_path: Path, out_dir: Path) -> None:
        return None

    def skip_reason(self, pt: MassPoint) -> str:
        """Why a requested point produced no row (for run_metadata.json)."""
        csv = self.vectors_path(pt)
        if not csv.exists() or csv.stat().st_size == 0:
            return "missing_or_empty_vectors"
        return "no_acceptance_or_nonpositive_ctau"


# ------------------------------------------------------------ the driver --

def run_point(spec: ModelSpec, pt: MassPoint, cfg: ScanConfig, mesh, *,
              rng=None, keep_arrays=False, return_mc=False) -> PointResult | None:
    """Acceptance + coupling scan for one mass point, or None when skipped."""
    pre = spec.pre_check(pt)
    if pre is not None:
        return PointResult(pre)

    csv_path = spec.vectors_path(pt)
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return None
    data = load_combined_csv(csv_path, pt.mass)
    n_events = len(data["weight"])
    if n_events == 0:
        return None

    hits, entry_d, exit_d = load_or_compute_geometry(
        spec.geometry_cache_path(pt, csv_path), data["eta"], data["phi"], mesh,
        force=cfg.force_geometry, source_mtime=csv_path.stat().st_mtime,
        batch_label=f"[{pt.tag}]")
    idx = np.where(hits & np.isfinite(entry_d) & np.isfinite(exit_d))[0]
    n_hits = int(hits.sum())
    if n_hits == 0:
        return PointResult(spec.empty_row(pt, n_events, 0, cfg))

    backend = spec.decay_backend(pt, csv_path)
    if backend is None:
        return None

    rng = rng if rng is not None else spec.rng_for(pt, cfg)
    idx, scan_weights, hit_estimator = spec.select_events(
        pt, idx, data["weight"][idx], cfg, rng)
    if len(idx) == 0:
        return PointResult(spec.empty_row(pt, n_events, n_hits, cfg, evaluated=True))
    scan_weights = spec.event_weights(pt, backend, scan_weights)

    entry_sel = entry_d[idx]
    exit_sel = exit_d[idx]
    path = exit_sel - entry_sel
    beta_gamma = data["beta_gamma"][idx]
    direction = directions_from_eta_phi(data["eta"][idx], data["phi"][idx])
    p_mag = beta_gamma * pt.mass
    energy = data["gamma"][idx] * pt.mass
    p4 = np.column_stack([energy, p_mag[:, None] * direction])

    grid = spec.grid.values()
    N_grid = np.zeros(len(grid))
    extra_scans = spec.extra_scans(pt, backend, cfg)
    extras = {e.tag: np.zeros(len(grid)) for e in extra_scans}
    d_parts, passed_parts, index_parts, mc_parts = [], [], [], []

    for sl, rng_c in spec.chunks(pt, len(idx), cfg, rng):
        d, passed, tmpl_idx, mc = backend.build(
            p4[sl], direction[sl], entry_sel[sl], exit_sel[sl], cfg.decay_samples, rng_c,
            return_mc=return_mc)
        d_parts.append(d)
        passed_parts.append(passed)
        index_parts.append(tmpl_idx)
        mc_parts.append(mc)
        if spec.scan_mode == "accumulate":
            _, N_part = scan_u2(d, passed, path[sl], scan_weights[sl], beta_gamma[sl],
                                backend.ctau_ref, L_INT_PB, grid)
            N_grid += N_part
            for extra in extra_scans:
                _, N_part = scan_u2(d, passed, path[sl], scan_weights[sl], beta_gamma[sl],
                                    backend.ctau_ref / extra.ctau_scale, L_INT_PB, grid,
                                    sample_w=extra.sample_weight(tmpl_idx))
                extras[extra.tag] += N_part

    d = np.concatenate(d_parts, axis=0)
    passed = np.concatenate(passed_parts, axis=0)
    template_index = (np.concatenate(index_parts, axis=0)
                      if all(p is not None for p in index_parts) else None)
    sample_w = backend.sample_weights(template_index) if template_index is not None else None
    if spec.scan_mode != "accumulate":
        _, N_grid = scan_u2(d, passed, path, scan_weights, beta_gamma,
                            backend.ctau_ref, L_INT_PB, grid, sample_w=sample_w)

    mc = None
    if return_mc and all(m is not None for m in mc_parts):
        mc = {k: np.concatenate([m[k] for m in mc_parts]) for k in mc_parts[0]}

    arrays = ScanArrays(grid=grid, N=N_grid, d=d, passed=passed, path=path,
                        weights=scan_weights, beta_gamma=beta_gamma,
                        ctau_ref=backend.ctau_ref, sample_w=sample_w, extras=extras,
                        n_hits_eval=len(idx), hit_estimator=hit_estimator,
                        n_samples=cfg.decay_samples, template_index=template_index, mc=mc)
    row = spec.finish(pt, spec.base_row(pt, n_events, n_hits), backend, arrays, cfg)
    return PointResult(row, arrays if (keep_arrays or return_mc) else None)


# ----------------------------------------------------------- the scan loop --

_WORKER = {}


def _worker_init(spec, cfg):
    _WORKER["mesh"] = get_mesh()
    _WORKER["spec"] = spec
    _WORKER["cfg"] = cfg


def _worker_point(pt):
    t0 = time.time()
    result = run_point(_WORKER["spec"], pt, _WORKER["cfg"], _WORKER["mesh"])
    return pt, result, time.time() - t0


def run_scan(spec: ModelSpec, points: list[MassPoint], cfg: ScanConfig, out_dir,
             *, workers: int = 1, resume: bool = False, verbose: bool = True) -> Path:
    """Scan ``points`` and write ``<out_dir>/sensitivity.csv``.

    With ``resume`` the rows already in that file are kept and their points
    skipped. Every completed point is checkpointed atomically, so an
    interrupted campaign loses at most the point in flight.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "sensitivity.csv"
    partial_csv = out_dir / "sensitivity.partial.csv"
    status_path = out_dir / "scan_status.json"

    rows: list[dict] = []
    done: set = set()
    if resume and out_csv.exists():
        previous = pd.read_csv(out_csv)
        rows = previous.to_dict("records")
        for r in rows:
            done.add((r.get("flavor") if spec.flavors else None, float(r["mass_GeV"])))
    todo = [pt for pt in points if pt.key not in done]
    if verbose:
        print(f"\n{spec.name}: {len(todo)} mass points to scan "
              f"({len(done)} resumed) with {workers} worker(s)...\n", flush=True)

    t_start = time.time()
    n_total = len(todo)

    def sort_key(r):
        return (r.get("flavor", "") or "", float(r["mass_GeV"]))

    def checkpoint(n_done):
        if rows:
            atomic_csv(pd.DataFrame(sorted(rows, key=sort_key)), partial_csv)
        elapsed = time.time() - t_start
        n_sens = sum(1 for r in rows if r.get("has_sensitivity"))
        atomic_json(status_path, {
            "ts": datetime.now().isoformat(), "done": n_done, "total": n_total,
            "n_sensitive": n_sens, "elapsed_s": round(elapsed, 1),
            "eta_s": round((elapsed / max(n_done, 1)) * (n_total - n_done), 1) if n_done else 0,
        }, sort_keys=False)

    def report(pt, result, elapsed, n_done):
        if not verbose:
            return
        if result is None:
            state = "skipped"
        elif result.row.get("exclusion_reason"):
            state = f"excluded ({result.row['exclusion_reason']})"
        elif result.row.get("has_sensitivity"):
            state = f"sensitive, peak_N={result.row.get('peak_N', float('nan')):.1f}"
        else:
            state = f"no sensitivity (peak_N={result.row.get('peak_N', float('nan')):.2f})"
        print(f"  [{n_done}/{n_total}] {pt.tag}: {state} ({elapsed:.1f}s)", flush=True)

    if workers <= 1:
        mesh = get_mesh()
        for i, pt in enumerate(todo, 1):
            t0 = time.time()
            result = run_point(spec, pt, cfg, mesh)
            if result is not None:
                rows.append(result.row)
            checkpoint(i)
            report(pt, result, time.time() - t0, i)
    else:
        with ProcessPoolExecutor(max_workers=workers, initializer=_worker_init,
                                 initargs=(spec, cfg)) as pool:
            futures = [pool.submit(_worker_point, pt) for pt in todo]
            for i, future in enumerate(as_completed(futures), 1):
                pt, result, elapsed = future.result()
                if result is not None:
                    rows.append(result.row)
                checkpoint(i)
                report(pt, result, elapsed, i)

    total_time = time.time() - t_start
    processed = {(r.get("flavor") if spec.flavors else None, float(r["mass_GeV"])) for r in rows}
    skipped = [{"flavor": pt.flavor, "mass_GeV": pt.mass, "reason": spec.skip_reason(pt)}
               for pt in points if pt.key not in processed]
    meta = {
        "model": spec.name,
        "timestamp": datetime.now().isoformat(),
        "flavors": list(spec.flavors) if spec.flavors else None,
        "masses_GeV": sorted({float(pt.mass) for pt in points}),
        "n_points_requested": len(points),
        "n_points_processed": len(rows),
        "n_points_skipped": len(skipped),
        "skipped": skipped,
        "workers": workers,
        "decay_samples": cfg.decay_samples,
        "thresholds": list(cfg.thresholds),
        "max_hit_events": cfg.max_hit_events,
        "event_chunk": cfg.event_chunk,
        "seed_salt": cfg.seed_salt,
        "seed_offset": cfg.seed_offset,
        "track_momentum_cut_GeV": P_CUT,
        "n_sensitive": sum(1 for r in rows if r.get("has_sensitivity")),
        "total_time_s": round(total_time, 1),
    }
    if skipped and verbose:
        by_reason: dict[str, int] = {}
        for s in skipped:
            by_reason[s["reason"]] = by_reason.get(s["reason"], 0) + 1
        print(f"\nWARNING: {len(skipped)}/{len(points)} requested points skipped:")
        for reason, n in sorted(by_reason.items()):
            print(f"    {n:4d}  {reason}")
    atomic_json(out_dir / "run_metadata.json", meta, sort_keys=False)
    if not rows:
        if verbose:
            print("\nNo results produced -- every requested point was skipped.")
        return out_csv
    atomic_csv(pd.DataFrame(sorted(rows, key=sort_key)), out_csv)
    if partial_csv.exists():
        partial_csv.unlink()
    if verbose:
        print(f"\nResults saved: {out_csv}")
        print(f"Sensitivity at {meta['n_sensitive']}/{len(rows)} mass points; "
              f"{total_time:.0f}s")
    spec.plot(out_csv, out_dir)
    return out_csv


# ------------------------------------------------------ threshold helpers --

def secondary_threshold_columns(result: dict, solve, secondary, fields,
                                rename=None) -> None:
    """Solve the same yield curve at each secondary threshold and add the
    ``<field>_N<T>`` columns to ``result``, in the model's field order."""
    for threshold in secondary:
        tag = f"N{threshold:g}"
        band = solve(threshold)
        for src in fields:
            dst = rename(src) if rename else src
            result[f"{dst}_{tag}"] = band[src]
