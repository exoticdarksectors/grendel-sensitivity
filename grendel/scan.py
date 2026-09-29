"""The per-mass sensitivity scan every benchmark model runs through."""
from __future__ import annotations

import functools
import hashlib
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
from .geometry.raycast import directions_from_eta_phi, get_mesh, load_or_compute_geometry, path_length
from .io.atomic import atomic_csv, atomic_json
from .io.vectors import format_mass_for_filename, load_combined_csv
from .reco.acceptance import P_CUT, scan_u2, signal_contribution_diagnostics


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
    thresholds: tuple[float, ...] = (3.0,)
    max_hit_events: int | None = None
    event_chunk: int | None = None
    seed_salt: str = ""
    seed_offset: int = 0
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
    """A co-varied yield curve accumulated alongside the nominal one (``scan_mode = "accumulate"``
    only): the lifetime divided by ``ctau_scale`` (a width scale) and each sample reweighted by
    ``sample_weight(template_index)``."""
    tag: str
    ctau_scale: float
    sample_weight: Callable[[np.ndarray], np.ndarray]


class DecayBackend(Protocol):
    """One mass point's decay model: how daughters are drawn, and the lifetime at the reference
    coupling (|U|^2 = 1, sin^2 theta = 1, 1/f = 1/f_ref) the scan reweights from."""
    name: str
    ctau_ref: float

    def build(self, p4, direction, entry_d, exit_d, n_samples, rng, return_mc=False):
        """-> (d, passed, template_index or None, mc or None); ``mc`` holds the full-length selection
        arrays (``acceptance.scatter_mc``) when asked for."""

    def sample_weights(self, template_index):
        """Per-sample weights (matrix-element reweighting) or None."""


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
    coupling_power: np.ndarray | None = None
    extras: dict[str, np.ndarray] = field(default_factory=dict)
    n_hits_eval: int = 0
    hit_estimator: str = "exact"
    n_samples: int = 0
    template_index: np.ndarray | None = None
    mc: dict | None = None

    def evaluate(self, u2: float) -> float:
        _, signal = scan_u2(self.d, self.passed, self.path, self.weights,
                            self.beta_gamma, self.ctau_ref, L_INT_PB,
                            np.asarray([u2]), sample_w=self.sample_w,
                            coupling_power=self.coupling_power)
        return signal[0]

    def diagnostics(self, u2: float) -> dict:
        return signal_contribution_diagnostics(
            self.d, self.passed, self.path, self.weights, self.beta_gamma,
            self.ctau_ref, u2, sample_w=self.sample_w,
            coupling_power=self.coupling_power)


@dataclass
class PointResult:
    row: dict
    arrays: ScanArrays | None = None


class ModelSpec:
    """What a benchmark model tells the driver."""

    name: str = ""
    grid: CouplingGrid = CouplingGrid(-12.0, -1.0, 200)
    scan_mode: str = "concatenate"
    decay_samples_default: int = 100
    flavors: tuple[str, ...] | None = None

    def points(self, masses, flavors=None) -> list[MassPoint]:
        if self.flavors is None:
            return [MassPoint(float(m)) for m in masses]
        return [MassPoint(float(m), f) for f in (flavors or self.flavors) for m in masses]

    def pre_check(self, pt: MassPoint) -> dict | None:
        """A result row that short-circuits the point, or None."""
        return None

    def vectors_path(self, pt: MassPoint) -> Path:
        raise NotImplementedError

    def input_files(self, pt: MassPoint) -> list[Path]:
        """The files the point's row is computed from; ``resume`` keeps a row only while they hash to
        what they did."""
        return [self.vectors_path(pt)]

    def geometry_cache_path(self, pt: MassPoint, vectors_path: Path) -> Path:
        raise NotImplementedError

    def decay_backend(self, pt: MassPoint, vectors_path: Path) -> DecayBackend | None:
        """None means the point is skipped (missing templates, closed decay)."""
        raise NotImplementedError

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

    def base_row(self, pt: MassPoint, n_events: int, n_hits: int) -> dict:
        raise NotImplementedError

    def empty_row(self, pt, n_events, n_hits, cfg, *, evaluated=False) -> dict:
        raise NotImplementedError

    def finish(self, pt, base: dict, backend, arrays: ScanArrays, cfg: ScanConfig) -> dict:
        raise NotImplementedError

    def skip_reason(self, pt: MassPoint) -> str:
        """Why a requested point produced no row (for run_metadata.json)."""
        csv = self.vectors_path(pt)
        if not csv.exists() or csv.stat().st_size == 0:
            return "missing_or_empty_vectors"
        return "no_acceptance_or_nonpositive_ctau"


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
    idx = np.where(hits & np.isfinite(entry_d[:, 0]) & np.isfinite(exit_d[:, 0]))[0]
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
    power = data["coupling_power"][idx]
    coupling_power = None if np.all(power == 1.0) else power

    entry_sel = entry_d[idx]
    exit_sel = exit_d[idx]
    path = path_length(entry_sel, exit_sel)
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
            power_sl = None if coupling_power is None else coupling_power[sl]
            sw = backend.sample_weights(tmpl_idx) if tmpl_idx is not None else None
            _, N_part = scan_u2(d, passed, path[sl], scan_weights[sl], beta_gamma[sl],
                                backend.ctau_ref, L_INT_PB, grid, sample_w=sw,
                                coupling_power=power_sl)
            N_grid += N_part
            for extra in extra_scans:
                ew = extra.sample_weight(tmpl_idx)
                _, N_part = scan_u2(d, passed, path[sl], scan_weights[sl], beta_gamma[sl],
                                    backend.ctau_ref / extra.ctau_scale, L_INT_PB, grid,
                                    sample_w=ew if sw is None else sw * ew,
                                    coupling_power=power_sl)
                extras[extra.tag] += N_part

    d = np.concatenate(d_parts, axis=0)
    passed = np.concatenate(passed_parts, axis=0)
    template_index = (np.concatenate(index_parts, axis=0)
                      if all(p is not None for p in index_parts) else None)
    sample_w = backend.sample_weights(template_index) if template_index is not None else None
    if spec.scan_mode != "accumulate":
        _, N_grid = scan_u2(d, passed, path, scan_weights, beta_gamma,
                            backend.ctau_ref, L_INT_PB, grid, sample_w=sample_w,
                            coupling_power=coupling_power)

    mc = None
    if return_mc and all(m is not None for m in mc_parts):
        mc = {k: np.concatenate([m[k] for m in mc_parts]) for k in mc_parts[0]}

    arrays = ScanArrays(grid=grid, N=N_grid, d=d, passed=passed, path=path,
                        weights=scan_weights, beta_gamma=beta_gamma,
                        ctau_ref=backend.ctau_ref, sample_w=sample_w,
                        coupling_power=coupling_power, extras=extras,
                        n_hits_eval=len(idx), hit_estimator=hit_estimator,
                        n_samples=cfg.decay_samples, template_index=template_index, mc=mc)
    row = spec.finish(pt, spec.base_row(pt, n_events, n_hits), backend, arrays, cfg)
    return PointResult(row, arrays if (keep_arrays or return_mc) else None)


_WORKER = {}


def _worker_init(spec, cfg):
    _WORKER["mesh"] = get_mesh()
    _WORKER["spec"] = spec
    _WORKER["cfg"] = cfg


def _worker_point(pt):
    t0 = time.time()
    digest = inputs_digest(_WORKER["spec"], pt)
    result = run_point(_WORKER["spec"], pt, _WORKER["cfg"], _WORKER["mesh"])
    return pt, result, digest, time.time() - t0


_CONFIG_FIELDS = ("decay_samples", "thresholds", "max_hit_events", "event_chunk",
                  "seed_salt", "seed_offset")
_PACKAGE = Path(__file__).resolve().parent


@functools.lru_cache(maxsize=1)
def code_digest() -> str:
    """SHA-256 over the package a scan runs: every file under ``grendel/`` (sources and package data)
    except ``band/``, which only post-processes scans."""
    h = hashlib.sha256()
    for path in sorted(_PACKAGE.rglob("*")):
        rel = path.relative_to(_PACKAGE)
        if (not path.is_file() or rel.parts[0] == "band" or path.suffix == ".pyc"
                or any(part.startswith(".") or part == "__pycache__" for part in rel.parts)):
            continue
        h.update(rel.as_posix().encode() + b"\0" + hashlib.sha256(path.read_bytes()).digest())
    return h.hexdigest()


def _file_sha256(path) -> str:
    h = hashlib.sha256()
    try:
        with open(path, "rb") as fh:
            for block in iter(lambda: fh.read(1 << 22), b""):
                h.update(block)
    except FileNotFoundError:
        return "missing"
    return h.hexdigest()


def inputs_digest(spec, pt: MassPoint) -> str:
    """SHA-256 over the contents of the files ``pt``'s row is computed from."""
    h = hashlib.sha256()
    for path in spec.input_files(pt):
        h.update(f"{Path(path).name}:{_file_sha256(path)}\n".encode())
    return h.hexdigest()


def _row_key(key: tuple) -> str:
    flavor, mass = key
    return f"{flavor or '-'}/{float(mass)!r}"


def _config_fingerprint(spec, cfg) -> dict:
    """What every row of a scan depends on besides its input files, for ``resume`` to compare: the
    model (name, grid, input directories, backend options), the ScanConfig and the package code."""
    model = {"name": spec.name, "grid": [spec.grid.log10_min, spec.grid.log10_max, spec.grid.n]}
    for key, value in sorted(vars(spec).items()):
        if key == "paths":
            model["vectors"], model["templates"] = str(value.vectors), str(value.templates)
        elif isinstance(value, Path):
            model[key] = str(value)
        elif value is None or isinstance(value, (str, int, float, bool)):
            model[key] = value
        else:
            model[key] = type(value).__name__
    scan = {k: getattr(cfg, k) for k in _CONFIG_FIELDS}
    scan["thresholds"] = [float(t) for t in scan["thresholds"]]
    return {"model": model, "scan": scan, "code": code_digest()}


def _resumable_rows(spec, config, out_csv, meta_path, partial_csv, status_path, verbose):
    """Rows of a finished (``sensitivity.csv``) and an interrupted (``sensitivity.partial.csv``) run
    that ``resume`` may keep, with their input digests: rows recorded under this configuration
    whose input files still hash to what they were computed from."""
    def recorded(path):
        try:
            return json.loads(Path(path).read_text())
        except (OSError, ValueError):
            return {}

    by_key: dict = {}
    digests: dict = {}
    for csv, meta_file in ((out_csv, meta_path), (partial_csv, status_path)):
        if not csv.exists():
            continue
        meta = recorded(meta_file)
        inputs = meta.get("inputs")
        if not isinstance(inputs, dict) or "config" not in meta:
            if verbose:
                print(f"  resume: {csv.name} records no input digests; its rows are recomputed")
            continue
        if meta["config"] != config:
            raise ValueError(f"resume: {csv} was produced with another configuration "
                             f"(recorded {meta['config']}, now {config}); rerun without --resume or "
                             "write to another --out")
        stale = 0
        for r in pd.read_csv(csv).to_dict("records"):
            pt = MassPoint(float(r["mass_GeV"]), r.get("flavor") if spec.flavors else None)
            k = _row_key(pt.key)
            if inputs.get(k) != inputs_digest(spec, pt):
                stale += 1
                continue
            by_key[pt.key] = r
            digests[k] = inputs[k]
        if stale and verbose:
            print(f"  resume: {stale} row(s) of {csv.name} came from other input files; recomputed")
    return list(by_key.values()), digests


class NoResultsError(RuntimeError):
    """A scan in which every requested point was skipped."""


def run_scan(spec: ModelSpec, points: list[MassPoint], cfg: ScanConfig, out_dir,
             *, workers: int = 1, resume: bool = False, verbose: bool = True) -> Path:
    """Scan ``points`` and write ``<out_dir>/sensitivity.csv``."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "sensitivity.csv"
    partial_csv = out_dir / "sensitivity.partial.csv"
    status_path = out_dir / "scan_status.json"

    config = _config_fingerprint(spec, cfg)
    rows: list[dict] = []
    digests: dict = {}
    done: set = set()
    if resume:
        rows, digests = _resumable_rows(spec, config, out_csv, out_dir / "run_metadata.json",
                                        partial_csv, status_path, verbose)
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
            "config": config,
            "inputs": digests,
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
            digest = inputs_digest(spec, pt)
            result = run_point(spec, pt, cfg, mesh)
            if result is not None:
                rows.append(result.row)
                digests[_row_key(pt.key)] = digest
            checkpoint(i)
            report(pt, result, time.time() - t0, i)
    else:
        with ProcessPoolExecutor(max_workers=workers, initializer=_worker_init,
                                 initargs=(spec, cfg)) as pool:
            futures = [pool.submit(_worker_point, pt) for pt in todo]
            for i, future in enumerate(as_completed(futures), 1):
                pt, result, digest, elapsed = future.result()
                if result is not None:
                    rows.append(result.row)
                    digests[_row_key(pt.key)] = digest
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
        "config": config,
        "inputs": digests,
    }
    by_reason: dict[str, int] = {}
    for s in skipped:
        by_reason[s["reason"]] = by_reason.get(s["reason"], 0) + 1
    if skipped and verbose:
        print(f"\nWARNING: {len(skipped)}/{len(points)} requested points skipped:")
        for reason, n in sorted(by_reason.items()):
            print(f"    {n:4d}  {reason}")
    atomic_json(out_dir / "run_metadata.json", meta, sort_keys=False)
    if not rows:
        moved = []
        for stale in (out_csv, partial_csv):
            if stale.exists():
                stale.replace(stale.with_name(stale.name + ".stale"))
                moved.append(stale.name)
        reasons = ", ".join(f"{n} {r}" for r, n in sorted(by_reason.items())) or "no points requested"
        raise NoResultsError(
            f"{spec.name}: no results -- every requested point was skipped ({reasons}); "
            f"see {out_dir / 'run_metadata.json'}"
            + (f"; the earlier {' and '.join(moved)} moved aside to *.stale" if moved else ""))
    atomic_csv(pd.DataFrame(sorted(rows, key=sort_key)), out_csv)
    if partial_csv.exists():
        partial_csv.unlink()
    if verbose:
        print(f"\nResults saved: {out_csv}")
        print(f"Sensitivity at {meta['n_sensitive']}/{len(rows)} mass points; "
              f"{total_time:.0f}s")
    return out_csv


def secondary_threshold_columns(result: dict, solve, secondary, fields,
                                rename=None) -> None:
    """Solve the same yield curve at each secondary threshold and add the ``<field>_N<T>`` columns to
    ``result``, in the model's field order."""
    for threshold in secondary:
        tag = f"N{threshold:g}"
        band = solve(threshold)
        for src in fields:
            dst = rename(src) if rename else src
            result[f"{dst}_{tag}"] = band[src]
