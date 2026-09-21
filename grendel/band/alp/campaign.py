"""Run and collect the independently generated BC10 uncertainty variations.

Every FONLL member and C_bs point receives a fresh production sample and its
own geometry/reconstruction/sensitivity tree; the gluon-surrogate variants
reuse the central production vectors (their production physics is
identical) but own independently simulated templates, geometry caches and
scans. Production and scan run as the package's own commands in a
per-variation environment.

Heavy artifacts live under the campaign directory (``--scratch-dir`` or
``GRENDEL_BC10_CAMPAIGN_DIR``); the runner is interruption-safe: mass CSVs
and scan checkpoints are atomic, and completion markers are written only
after the output set validates and hashes.

    python -m grendel.band.alp.campaign run --grid-dir G --scratch-dir S \\
        --central-template-dir T --axes fonll cbs numerical_control
    python -m grendel.band.alp.campaign collect --grid-dir G --scratch-dir S \\
        --central-curve RESULTS/bc10/sensitivity.csv --out-dir RESULTS/bc10/band
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from ...io.atomic import atomic_csv, atomic_json, sha256_file
from ...io.paths import repo_root
from ...io.vectors import format_mass_for_filename
from ...models.alp import model
from ...models.alp.mass_grid import ALP_MASS_GRID
from ..campaign_store import git_state, remove_generated_tree, tree_usage
from .combine import combine_band
from .variations import all_variations


def sha256_tree(paths) -> str:
    """Hash file names and contents in stable order."""
    digest = hashlib.sha256()
    for path in sorted((Path(p) for p in paths), key=lambda p: str(p)):
        digest.update(path.name.encode())
        digest.update(b"\0")
        digest.update(sha256_file(path).encode())
        digest.update(b"\n")
    return digest.hexdigest()


# ------------------------------------------------------------- mass grid --

def load_mass_grid(path: Path | None) -> dict:
    """The exact mass list used by every campaign stage, pinned by hash."""
    if path is None:
        masses = [float(v) for v in ALP_MASS_GRID]
        source, source_path, source_sha256 = "grendel.models.alp.mass_grid.ALP_MASS_GRID", None, None
    else:
        path = Path(path).expanduser().resolve()
        frame = pd.read_csv(path)
        if "mass_GeV" not in frame:
            raise ValueError(f"mass-grid file has no mass_GeV column: {path}")
        masses = frame["mass_GeV"].astype(float).tolist()
        source, source_path, source_sha256 = "csv", str(path), sha256_file(path)
    values = np.asarray(masses, dtype=float)
    if len(values) == 0 or np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("mass grid must contain finite positive values")
    if np.any(np.diff(values) <= 0.0):
        raise ValueError("mass grid must be strictly increasing with no duplicates")
    canonical = json.dumps(masses, separators=(",", ":")).encode()
    return {"source": source, "source_path": source_path, "source_sha256": source_sha256,
            "canonical_sha256": hashlib.sha256(canonical).hexdigest(),
            "n_masses": len(masses), "masses_GeV": masses}


def supported_masses(masses) -> list[float]:
    return [float(m) for m in masses if model.excluded_light_meson_resonance(m) is None]


# ------------------------------------------------------------ validation --

def _expected_vector_paths(directory: Path, masses):
    return [directory / f"mA_{format_mass_for_filename(m)}.csv" for m in supported_masses(masses)]


def _expected_template_paths(directory: Path, masses):
    return [directory / f"templates_{format_mass_for_filename(m)}.npz" for m in supported_masses(masses)]


def validate_vectors(directory: Path, masses, hash_outputs=True) -> dict:
    paths = _expected_vector_paths(directory, masses)
    missing = [str(p) for p in paths if not p.exists() or p.stat().st_size == 0]
    if missing:
        raise RuntimeError(f"production tree is incomplete ({len(missing)} missing/empty); first: {missing[0]}")
    temporary = list(directory.glob("*.tmp"))
    if temporary:
        raise RuntimeError(f"production tree has partial temporary file: {temporary[0]}")
    return {"n_vector_files": len(paths),
            "total_vector_bytes": sum(p.stat().st_size for p in paths),
            "vector_tree_sha256": sha256_tree(paths) if hash_outputs else None}


def validate_sensitivity_csv(path: Path, masses) -> dict:
    if not path.exists():
        raise RuntimeError(f"missing sensitivity CSV: {path}")
    frame = pd.read_csv(path)
    expected = [float(v) for v in masses]
    observed = sorted(frame["mass_GeV"].astype(float))
    if len(frame) != len(expected) or observed != sorted(expected):
        raise RuntimeError(f"sensitivity checkpoint has {len(frame)} rows; expected {len(expected)}")
    if path.with_suffix(path.suffix + ".tmp").exists() or (path.parent / "sensitivity.partial.csv").exists():
        raise RuntimeError(f"partial sensitivity checkpoint remains next to {path}")
    return {"sensitivity_csv_sha256": sha256_file(path), "n_sensitivity_rows": len(frame),
            "n_sensitive_rows": int(frame["has_sensitivity"].sum())}


def validate_geometry(geometry_dir: Path, masses) -> dict:
    geometry = sorted(geometry_dir.glob("geom_*.npz"))
    expected = len(supported_masses(masses))
    if len(geometry) != expected:
        raise RuntimeError(f"geometry cache has {len(geometry)} files; expected {expected}")
    return {"n_geometry_files": len(geometry),
            "total_geometry_bytes": sum(p.stat().st_size for p in geometry),
            "geometry_tree_sha256": sha256_tree(geometry)}


def template_provenance(directory: Path, masses, expected_decay_model=None) -> dict:
    paths = _expected_template_paths(directory, masses)
    missing = [p for p in paths if not p.exists()]
    if missing:
        raise RuntimeError(f"template directory {directory} is missing {len(missing)} campaign files; "
                           f"first: {missing[0]}")
    provenance = {"path": str(directory.resolve()), "n_files": len(paths),
                  "tree_sha256": sha256_tree(paths),
                  "total_bytes": sum(p.stat().st_size for p in paths)}
    if expected_decay_model is not None:
        observed, counts = set(), set()
        for p in paths:
            with np.load(p) as bundle:
                if "decay_model" not in bundle.files:
                    raise RuntimeError(f"{p} has no pinned decay_model metadata")
                observed.add(str(bundle["decay_model"]))
                counts.add(int(bundle["n_templates"]))
        if observed != {expected_decay_model}:
            raise RuntimeError(f"template decay models {observed} do not match {expected_decay_model}")
        provenance.update({"decay_model": expected_decay_model, "n_templates_per_mass": sorted(counts)})
    return provenance


def _run_logged(command, env, log_path: Path, label: str) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"[{label}] {' '.join(command)}", flush=True)
    with log_path.open("a") as log:
        log.write(f"\n=== {time.strftime('%Y-%m-%dT%H:%M:%S%z')} ===\n$ " + " ".join(command) + "\n")
        log.flush()
        result = subprocess.run(command, cwd=repo_root(), env=env, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"{label} failed with exit {result.returncode}; see {log_path}")


def _marker_matches(path: Path, variation, n_pool, code, mass_grid) -> bool:
    if not path.exists():
        return False
    payload = json.loads(path.read_text())
    return (payload.get("variation", {}).get("name") == variation["name"]
            and payload.get("variation", {}).get("grid_sha256") == variation["grid_sha256"]
            and payload.get("n_pool") == n_pool
            and payload.get("code", {}) == code
            and payload.get("mass_grid", {}).get("canonical_sha256") == mass_grid["canonical_sha256"])


def _validate_completed_run(run_dir: Path, marker: dict) -> None:
    masses = marker["mass_grid"]["masses_GeV"]
    current = validate_sensitivity_csv(Path(marker["sensitivity_csv"]), masses)
    if current["sensitivity_csv_sha256"] != marker["sensitivity_csv_sha256"]:
        raise RuntimeError(f"completed {marker['variation']['name']} sensitivity checksum changed")
    state = marker.get("storage_state", "full")
    if state == "full":
        geometry = validate_geometry(run_dir / "analysis" / "geometry_cache", masses)
        if geometry["geometry_tree_sha256"] != marker["geometry_tree_sha256"]:
            raise RuntimeError(f"completed {marker['variation']['name']} geometry checksum changed")
        if marker["variation"]["production_mode"] == "fresh_600k":
            vectors = validate_vectors(run_dir / "llp_4vectors", masses)
            production = json.loads(Path(marker["production_marker"]).read_text())
            if vectors["vector_tree_sha256"] != production["vector_tree_sha256"]:
                raise RuntimeError(f"completed {marker['variation']['name']} vector checksum changed")
    elif state == "compacted":
        for path in (run_dir / "llp_4vectors", run_dir / "analysis" / "geometry_cache"):
            if path.exists():
                raise RuntimeError(f"compacted run unexpectedly retains generated tree: {path}")
    elif state != "compacting":
        raise RuntimeError(f"unknown completion storage state: {state}")


def _compact_completed_run(completion: Path, marker: dict, run_dir: Path, keep_intermediates: bool) -> dict:
    variation = marker["variation"]
    state = marker.get("storage_state", "full")
    if variation["name"] == "central" or (keep_intermediates and state == "full") or state == "compacted":
        return marker
    if state not in {"full", "compacting"}:
        raise RuntimeError(f"cannot compact completion in storage state {state}")
    vector_dir = run_dir / "llp_4vectors"
    geometry_dir = run_dir / "analysis" / "geometry_cache"
    targets = [geometry_dir]
    if variation["production_mode"] == "fresh_600k":
        targets.insert(0, vector_dir)
    if state == "full":
        marker = dict(marker)
        marker["storage_state"] = "compacting"
        marker["compaction"] = {"started_unix": time.time(),
                                "targets": [tree_usage(p) for p in targets if p.exists()],
                                "retained": [marker["sensitivity_csv"], marker["production_marker"],
                                             str(run_dir / "production.log"), str(run_dir / "sensitivity.log")]}
        atomic_json(completion, marker)
    for target in targets:
        remove_generated_tree(target, run_dir)
    marker = dict(marker)
    marker["storage_state"] = "compacted"
    marker["compaction"] = dict(marker["compaction"])
    marker["compaction"]["completed_unix"] = time.time()
    marker["compaction"]["reclaimed_bytes"] = sum(t["bytes"] for t in marker["compaction"]["targets"])
    atomic_json(completion, marker)
    print(f"[{variation['name']}] compacted {marker['compaction']['reclaimed_bytes'] / 2**30:.1f} GiB", flush=True)
    return marker


def _environment(run_dir, vector_dir, analysis_dir, template_dir, variation) -> dict:
    env = os.environ.copy()
    env.update({
        "GRENDEL_WORK_DIR": str(run_dir),
        "GRENDEL_BC10_VECTORS_DIR": str(vector_dir),
        "GRENDEL_BC10_ANALYSIS_DIR": str(analysis_dir),
        "GRENDEL_BC10_GEOMETRY_DIR": str(analysis_dir / "geometry_cache"),
        "GRENDEL_BC10_TEMPLATES_DIR": str(template_dir),
        "GRENDEL_FONLL_BOTTOM_GRID": variation["grid_path"],
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
    })
    return env


def _validate_central_production_reference(marker_path: Path, vector_dir: Path, variation, n_pool,
                                           hash_vectors: bool, mass_grid: dict):
    marker = json.loads(marker_path.read_text())
    central = marker.get("variation", {})
    if central.get("name") != "central":
        raise RuntimeError(f"central production marker is mislabeled: {marker_path}")
    if central.get("grid_sha256") != variation["grid_sha256"]:
        raise RuntimeError("reused central production has a different FONLL grid")
    if marker.get("n_pool") != n_pool:
        raise RuntimeError("reused central production has a different pool size")
    if float(central.get("cbs_amplitude_scale", 1.0)) != 1.0:
        raise RuntimeError("reused central production has non-central C_bs")
    if marker.get("mass_grid", {}).get("canonical_sha256") != mass_grid["canonical_sha256"]:
        raise RuntimeError("reused central production has a different mass grid")
    info = validate_vectors(vector_dir, mass_grid["masses_GeV"], hash_outputs=hash_vectors)
    if hash_vectors and info["vector_tree_sha256"] != marker.get("vector_tree_sha256"):
        raise RuntimeError("reused central production vector checksum changed")
    info["vector_tree_sha256"] = marker.get("vector_tree_sha256")
    return marker, info


# --------------------------------------------------------------- run side --

def run_variation(variation, args, code, template_cache) -> None:
    mass_grid = args.mass_grid
    masses = mass_grid["masses_GeV"]
    scratch = Path(args.scratch_dir).expanduser().resolve()
    run_dir = scratch / "runs" / variation["name"]
    run_dir.mkdir(parents=True, exist_ok=True)
    analysis_dir = run_dir / "analysis"
    central_dir = scratch / "runs" / "central"
    if variation["template_variant"] == "central":
        template_dir = Path(args.central_template_dir).expanduser().resolve()
    else:
        template_dir = Path(args.gluon_template_root).expanduser().resolve() / variation["template_variant"]
    key = str(template_dir)
    if key not in template_cache:
        template_cache[key] = template_provenance(template_dir, masses)

    completion = run_dir / "variation.complete.json"
    if _marker_matches(completion, variation, args.n_pool, code, mass_grid):
        marker = json.loads(completion.read_text())
        _validate_completed_run(run_dir, marker)
        _compact_completed_run(completion, marker, run_dir, args.keep_intermediates)
        print(f"[{variation['name']}] variation already complete", flush=True)
        return
    if completion.exists():
        raise RuntimeError(f"stale completion marker for {variation['name']}; use a new scratch directory")

    if variation["production_mode"] == "fresh_600k":
        vector_dir = run_dir / "llp_4vectors"
        production_marker = run_dir / "production.complete.json"
        if _marker_matches(production_marker, variation, args.n_pool, code, mass_grid):
            info = validate_vectors(vector_dir, masses, hash_outputs=True)
            if info["vector_tree_sha256"] != json.loads(production_marker.read_text())["vector_tree_sha256"]:
                raise RuntimeError(f"completed production checksum changed for {variation['name']}")
            print(f"[{variation['name']}] production already complete", flush=True)
        else:
            if production_marker.exists():
                raise RuntimeError(f"stale production marker for {variation['name']}; use a new scratch "
                                   "directory or remove that generated run explicitly")
            env = _environment(run_dir, vector_dir, analysis_dir, template_dir, variation)
            command = [sys.executable, "-u", "-m", "grendel.models.alp.production",
                       "--n-pool", str(args.n_pool), "--seed", str(variation["production_seed"]),
                       "--cbs-amplitude-scale", str(variation["cbs_amplitude_scale"]),
                       "--out-dir", str(vector_dir), "--mass", *(str(m) for m in masses), "--resume"]
            started = time.time()
            log = run_dir / "production.log"
            _run_logged(command, env, log, f"{variation['name']}:production")
            info = validate_vectors(vector_dir, masses)
            atomic_json(production_marker, {"variation": variation, "n_pool": args.n_pool,
                                            "mass_grid": mass_grid, "code": code, "command": command,
                                            "log": str(log), "log_sha256": sha256_file(log),
                                            "started_unix": started, "completed_unix": time.time(), **info})
    else:
        central_marker = central_dir / "production.complete.json"
        if not central_marker.exists():
            raise RuntimeError(f"{variation['name']} requires completed central production: {central_marker}")
        vector_dir = central_dir / "llp_4vectors"
        _, info = _validate_central_production_reference(central_marker, vector_dir, variation, args.n_pool,
                                                         hash_vectors=False, mass_grid=mass_grid)
        atomic_json(run_dir / "production.reference.json", {
            "variation": variation, "central_production_marker": str(central_marker),
            "central_production_marker_sha256": sha256_file(central_marker),
            "vector_dir": str(vector_dir), "mass_grid": mass_grid, **info})

    env = _environment(run_dir, vector_dir, analysis_dir, template_dir, variation)
    output = analysis_dir / "sensitivity.csv"
    command = [sys.executable, "-u", "-m", "grendel.models.alp.scan", "--out", str(analysis_dir),
               "--vectors-dir", str(vector_dir), "--templates-dir", str(template_dir),
               "--reco-seed-offset", str(variation["reco_seed_offset"]),
               "--mass", *(str(m) for m in masses), "--resume", "--thresholds", *map(str, args.thresholds)]
    started = time.time()
    log = run_dir / "sensitivity.log"
    _run_logged(command, env, log, f"{variation['name']}:sensitivity")
    info = {**validate_sensitivity_csv(output, masses), **validate_geometry(analysis_dir / "geometry_cache", masses)}
    production_marker = (run_dir / "production.complete.json" if variation["production_mode"] == "fresh_600k"
                         else run_dir / "production.reference.json")
    atomic_json(completion, {
        "variation": variation, "n_pool": args.n_pool, "mass_grid": mass_grid, "code": code,
        "python": {"executable": sys.executable, "version": platform.python_version()},
        "template": template_cache[key],
        "production_marker": str(production_marker), "production_marker_sha256": sha256_file(production_marker),
        "sensitivity_csv": str(output), "command": command,
        "sensitivity_log": str(log), "sensitivity_log_sha256": sha256_file(log),
        "storage_state": "full",
        "regeneration": {"grid_path": variation["grid_path"], "grid_sha256": variation["grid_sha256"],
                         "production_seed": variation["production_seed"],
                         "reco_seed_offset": variation["reco_seed_offset"], "n_pool": args.n_pool,
                         "cbs_amplitude_scale": variation["cbs_amplitude_scale"],
                         "production_mode": variation["production_mode"],
                         "template_path": template_cache[key]["path"],
                         "template_tree_sha256": template_cache[key]["tree_sha256"],
                         "sensitivity_command": command},
        "started_unix": started, "completed_unix": time.time(), **info})
    _compact_completed_run(completion, json.loads(completion.read_text()), run_dir, args.keep_intermediates)
    print(f"[{variation['name']}] COMPLETE", flush=True)


# ----------------------------------------------------------- collect side --

def _scratch_relative(path: Path, scratch: Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(scratch.resolve()))
    except ValueError as exc:
        raise ValueError(f"campaign artifact is outside the campaign directory: {path}") from exc


def load_completed(scratch: Path, variations):
    """The raw variation table and a portable registry of the completed runs."""
    scratch = Path(scratch).expanduser().resolve()
    rows, registry = [], []
    for variation in variations:
        marker_path = scratch / "runs" / variation["name"] / "variation.complete.json"
        if not marker_path.exists():
            raise FileNotFoundError(f"variation {variation['name']} is incomplete: {marker_path}")
        marker = json.loads(marker_path.read_text())
        if marker["variation"] != variation:
            raise ValueError(f"variation definition mismatch in {marker_path}")
        sensitivity = Path(marker["sensitivity_csv"])
        actual_sha = sha256_file(sensitivity)
        if actual_sha != marker["sensitivity_csv_sha256"]:
            raise ValueError(f"sensitivity checksum mismatch: {sensitivity}")
        frame = pd.read_csv(sensitivity)
        expected_rows = int(marker["mass_grid"]["n_masses"])
        if len(frame) != expected_rows:
            raise ValueError(f"{sensitivity} has {len(frame)} rows, expected {expected_rows}")
        frame.insert(0, "axis", variation["axis"])
        frame.insert(0, "variation", variation["name"])
        rows.append(frame)
        portable = dict(variation)
        portable["grid_file"] = Path(portable.pop("grid_path")).name
        registry.append({
            "variation": portable,
            "completion_marker": _scratch_relative(marker_path, scratch),
            "completion_marker_sha256": sha256_file(marker_path),
            "sensitivity_csv": _scratch_relative(sensitivity, scratch),
            "sensitivity_csv_sha256": actual_sha,
            "production_marker": _scratch_relative(Path(marker["production_marker"]), scratch),
            "production_marker_sha256": marker["production_marker_sha256"],
            "template": {k: v for k, v in marker["template"].items() if k != "path"}
            | {"path_role": variation["template_variant"]},
            "code": marker["code"],
            "storage_state": marker.get("storage_state", "full"),
            "pre_compaction_hashes": {
                "vector_tree_sha256": json.loads(Path(marker["production_marker"]).read_text()).get("vector_tree_sha256"),
                "geometry_tree_sha256": marker["geometry_tree_sha256"]},
        })
    return pd.concat(rows, ignore_index=True), registry


def _topology(band: pd.DataFrame, stem: str) -> dict:
    """Where a variation family changes the island topology, and which member."""
    differs = band[f"{stem}_differs"].astype(bool)
    masses = [float(v) for v in band.loc[differs, "mass_GeV"]]
    names = [str(v) for v in band.loc[differs, f"{stem}_difference_variations"]]
    return {"n_differences": int(differs.sum()), "difference_masses_GeV": masses,
            "difference_variations_by_mass": {f"{m:g}": n for m, n in zip(masses, names)}}


def headline(band: pd.DataFrame) -> dict:
    sensitive = band[band["has_sensitivity"]].copy()
    out = {"n_mass_points": len(band), "n_central_sensitive": int(band["has_sensitivity"].sum()),
           "n_any_variation_sensitive": int(band["any_variation_sensitive"].sum()),
           "numerical_control_topology": _topology(band, "numerical_control_topology"),
           "physical_variation_topology": {
               **_topology(band, "halo_topology"),
               "restored_masses_GeV": [float(v) for v in band.loc[band["halo_restores_sensitivity"], "mass_GeV"]],
               "removed_masses_GeV": [float(v) for v in band.loc[band["halo_removes_sensitivity"], "mass_GeV"]]}}
    for boundary in ("invf_min", "invf_max"):
        central = sensitive[f"{boundary}_central"]
        for direction, edge in (("up", "hi"), ("dn", "lo")):
            envelope = sensitive[f"{boundary}_envelope_{edge}"]
            values = (np.log10(envelope / central) if direction == "up"
                      else np.log10(central / envelope)).replace([np.inf, -np.inf], np.nan)
            if values.notna().any():
                index = values.idxmax()
                out[f"max_{boundary}_envelope_{direction}_dex"] = float(values.loc[index])
                out[f"max_{boundary}_envelope_{direction}_dex_mass_GeV"] = float(band.loc[index, "mass_GeV"])
        repeat = sensitive[f"{boundary}_repeat_max_abs_dex"].replace([np.inf, -np.inf], np.nan)
        flagged = sensitive[sensitive[f"{boundary}_repeat_not_subdominant"]]
        out[f"{boundary}_numerical_control"] = {
            "median_repeat_max_abs_dex": float(repeat.median()) if repeat.notna().any() else None,
            "max_repeat_abs_dex": float(repeat.max()) if repeat.notna().any() else None,
            "n_not_subdominant": len(flagged),
            "not_subdominant_masses_GeV": [float(v) for v in flagged["mass_GeV"]]}
    return out


def collect(scratch: Path, variations, central_curve: Path, out_dir: Path) -> None:
    raw, registry = load_completed(scratch, variations)
    band = combine_band(raw, central_curve)
    out_dir = Path(out_dir)
    raw_path = out_dir / "bc10_uncertainty_variations.csv"
    band_path = out_dir / "bc10_single_source_variation_envelope.csv"
    atomic_csv(raw, raw_path)
    atomic_csv(band, band_path)
    codes = {json.dumps(item["code"], sort_keys=True) for item in registry}
    if len(codes) != 1:
        raise ValueError("variations were produced from different code states")
    atomic_json(out_dir / "UNCERTAINTY_MANIFEST.json", {
        "artifact": "GRENDEL BC10 single-source variation envelope",
        "generated_unix": time.time(),
        "method": {
            "production": "independent production, geometry, reconstruction and scan for the central "
                          "plus 6 scale, 100 NNPDF replica and 2 bottom-mass FONLL grids; no reweighting",
            "decay_gg": "independently simulated full-branching Pythia u/d/s surrogate templates, each "
                        "propagated through its own geometry/reconstruction run",
            "cbs": "+/-20% C_bs amplitude, implemented as fresh production rates scaled by 0.8^2 and 1.2^2",
            "numerical_control": "two same-physics central repeats with fresh production pools and distinct "
                                 "reconstruction seeds; excluded from the physical envelope",
            "combination": "pointwise one-source-at-a-time intervals in log10(1/f): named scale/mb/gg/C_bs "
                           "extrema and NNPDF replica 16th/84th percentiles; display envelope is their "
                           "outermost boundary; numerical repeats excluded; no quadrature, no "
                           "confidence-interval claim",
            "rebase": "one-source log10(1/f) shifts applied to the canonical central contour; numerical "
                      "repeats retain their directly simulated absolute values",
            "label": "single_source_variation_envelope"},
        "variation_counts": {axis: int(sum(v["axis"] == axis for v in variations))
                             for axis in ("central", "scale", "pdf", "mb", "decay_gg", "cbs", "numerical_control")},
        "inputs": {"canonical_central_curve": {"file": Path(central_curve).name,
                                               "sha256": sha256_file(central_curve)}},
        "outputs": {raw_path.name: sha256_file(raw_path), band_path.name: sha256_file(band_path)},
        "headline": headline(band),
        "code": json.loads(next(iter(codes))),
        "registry": registry})
    print(f"wrote {raw_path} ({len(raw)} rows), {band_path} ({len(band)} rows) and the manifest")


# --------------------------------------------------------------------- CLI --

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)

    def common(p):
        p.add_argument("--grid-dir", type=Path, default=os.environ.get("GRENDEL_FONLL_GRID_DIR"),
                       required=os.environ.get("GRENDEL_FONLL_GRID_DIR") is None)
        p.add_argument("--scratch-dir", type=Path, default=os.environ.get("GRENDEL_BC10_CAMPAIGN_DIR"),
                       required=os.environ.get("GRENDEL_BC10_CAMPAIGN_DIR") is None)
        p.add_argument("--mass-grid-file", type=Path, default=None,
                       help="CSV with a strictly increasing mass_GeV column (default: the BC10 grid)")

    run = sub.add_parser("run", help="run the selected variations")
    common(run)
    run.add_argument("--central-template-dir", type=Path, required=True)
    run.add_argument("--gluon-template-root", type=Path, default=None,
                     help="root containing gg_u/gg_d/gg_s template directories (decay_gg axis)")
    run.add_argument("--axes", nargs="+", choices=("fonll", "decay_gg", "cbs", "numerical_control"),
                     default=["fonll"])
    run.add_argument("--only", nargs="+", default=None)
    run.add_argument("--n-pool", type=int, default=600_000)
    run.add_argument("--worker-index", type=int, default=0)
    run.add_argument("--worker-count", type=int, default=1)
    run.add_argument("--max-variations", type=int, default=None)
    run.add_argument("--thresholds", nargs="+", type=float, default=[3.0])
    run.add_argument("--keep-intermediates", action="store_true",
                     help="retain non-central vectors and geometry caches after completion")
    coll = sub.add_parser("collect", help="validate and combine the completed campaign")
    common(coll)
    coll.add_argument("--central-curve", type=Path, required=True,
                      help="canonical BC10 sensitivity.csv used to rebase the campaign shifts")
    coll.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)

    try:
        args.mass_grid = load_mass_grid(args.mass_grid_file)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        ap.error(str(exc))
    variations = all_variations(args.grid_dir)
    if args.command == "collect":
        args.out_dir.mkdir(parents=True, exist_ok=True)
        collect(args.scratch_dir, variations, args.central_curve, args.out_dir)
        return 0
    if not 0 <= args.worker_index < args.worker_count:
        ap.error("worker-index must satisfy 0 <= index < worker-count")
    if "decay_gg" in args.axes and args.gluon_template_root is None:
        ap.error("--gluon-template-root is required for the decay_gg axis")
    Path(args.scratch_dir).expanduser().resolve().mkdir(parents=True, exist_ok=True)
    selected = [v for v in variations if v["campaign_axis"] in args.axes]
    if args.only:
        wanted = set(args.only)
        selected = [v for v in selected if v["name"] in wanted]
        missing = wanted - {v["name"] for v in selected}
        if missing:
            ap.error(f"unknown/ineligible variations: {sorted(missing)}")
    selected = [v for i, v in enumerate(selected) if i % args.worker_count == args.worker_index]
    if args.max_variations is not None:
        selected = selected[:args.max_variations]
    code = git_state()
    print(f"BC10 uncertainty worker {args.worker_index}/{args.worker_count}: {len(selected)} variations, "
          f"{args.mass_grid['n_masses']} masses under {args.scratch_dir}", flush=True)
    template_cache = {}
    for variation in selected:
        run_variation(variation, args, code, template_cache)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
