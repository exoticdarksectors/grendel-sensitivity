"""Run and collect the independent full-statistics BC4 uncertainty campaign.

Every variation receives its own freshly sampled FONLL parent pool, scalar
four-vectors, ray casts, decay/reconstruction Monte Carlo and coupling scan.
No events, weights, geometry or reconstruction outcomes are shared between
variations (``variations`` lists them).

Heavy products live under the campaign directory (``--scratch-dir`` or
``GRENDEL_BC4_CAMPAIGN_DIR``)::

    runs/<variation>/llp_4vectors/*.csv
    runs/<variation>/geometry_cache/*.npz
    runs/<variation>/results/*.json

Each four-vector file, mass result and aggregate curve is committed
atomically. After a variation passes full checksum validation its raw
vectors and geometry are reclaimed; compact results, tree hashes, seeds,
provenance and logs remain. ``run`` resumes incomplete variations and
recognises compacted ones; ``collect`` requires every variation to be
compacted and writes the raw table, the display envelope
(``combine.combine_band``) and the provenance manifest.

    python -m grendel.band.scalar.campaign status  --grid-dir G --scratch-dir S
    python -m grendel.band.scalar.campaign run     --grid-dir G --scratch-dir S --workers 2
    python -m grendel.band.scalar.campaign collect --grid-dir G --scratch-dir S \\
        --central-curve RESULTS/bc4/sensitivity.csv --out-dir RESULTS/bc4/band
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from functools import lru_cache
from importlib.metadata import PackageNotFoundError, version as package_version
from pathlib import Path

import numpy as np
import pandas as pd

from ...io.atomic import atomic_csv, atomic_json, sha256_file, tree_hash
from ...io.paths import ModelPaths
from ...models.scalar import production
from ...models.scalar.spec import ScalarSpec
from ...production.fonll.fonll_parser import get_sigma_total
from ...production.fonll.meson_sampler import sample_meson_4vectors
from ...scan import MassPoint, ScanConfig, run_point
from ..campaign_store import (append_log, code_hashes, config_sha256, file_record, git_head,
                              jsonable, read_json, software, stable_seed, stage_record, utc_now,
                              validate_recorded_code_state, variation_lock)
from .combine import combine_band
from .variations import NUMERICAL_CONTROL_VARIATIONS, discover_variations, select_variations

LOG_S2T_MIN, LOG_S2T_MAX, N_S2T = ScalarSpec.grid.log10_min, ScalarSpec.grid.log10_max, ScalarSpec.grid.n
GRID_ENV = "GRENDEL_FONLL_BOTTOM_GRID"
BOUNDARIES = (("u2_min", "u2_min_open"), ("u2_max", "u2_max_open"))

# Source files whose hashes are recorded with every run; a collector proves
# they match the producer commit before trusting retained results.
CODE_INPUTS = (
    "grendel/models/scalar/model.py",
    "grendel/models/scalar/production.py",
    "grendel/models/scalar/acceptance.py",
    "grendel/models/scalar/spec.py",
    "grendel/band/scalar/campaign.py",
    "grendel/scan.py",
    "grendel/reco/acceptance.py",
    "grendel/reco/exclusion.py",
    "grendel/geometry/raycast.py",
    "grendel/geometry/reco_common.py",
    "grendel/geometry/grendel_geometry.py",
)
CONFIG_KEYS = ("schema_version", "variation", "axis", "width_scheme", "grid_file", "grid_sha256",
               "n_parent_pool", "n_scalar_events_expected", "n_decay_samples_per_hit", "base_seed",
               "parent_pool_seed", "mass_grid_GeV", "coupling_grid", "ray_backend", "code_sha256")

RAY_BACKEND_VALIDATION = {
    "sample": "25000 scalar rays from the BC4 production pipeline",
    "reference": "trimesh.ray.ray_triangle.RayMeshIntersector",
    "hit_mask_mismatches": 0,
    "max_entry_distance_difference_m": 1.42e-14,
    "max_exit_distance_difference_m": 2.13e-14,
}


# ------------------------------------------------------------ provenance --

@lru_cache(maxsize=None)
def ray_backend_provenance() -> dict:
    """The pinned Embree ray backend the full campaign requires."""
    import trimesh
    from ...geometry.raycast import get_mesh
    mesh = get_mesh()
    backend = f"{type(mesh.ray).__module__}.{type(mesh.ray).__name__}"
    try:
        embreex_version = package_version("embreex")
    except PackageNotFoundError:
        embreex_version = None
    if embreex_version != "4.4.0" or "ray_pyembree" not in backend:
        raise RuntimeError("the BC4 full campaign requires the pinned Embree ray backend; "
                           "install embreex==4.4.0 in the active environment")
    return {"implementation": backend, "trimesh_version": trimesh.__version__,
            "embreex_version": embreex_version,
            "validation_against_triangle_backend": RAY_BACKEND_VALIDATION}


@lru_cache(maxsize=None)
def templates_provenance(templates_dir: str) -> dict:
    """Identity of a cached decay-template set: its generator manifest and a
    digest over every bundle, so a template-driven variation's retained
    config pins the sampled representation of the decay model."""
    root = Path(templates_dir)
    manifest_path = root / "MANIFEST.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"no template manifest under {root}")
    manifest = json.loads(manifest_path.read_text())
    bundles = sorted(root.glob("templates_*.npz"))
    if not bundles:
        raise FileNotFoundError(f"no template bundles under {root}")
    entries = [{"path": b.name, "bytes": b.stat().st_size, "sha256": sha256_file(b)} for b in bundles]
    return {"directory": root.name, "manifest_sha256": sha256_file(manifest_path),
            "bundle_tree_sha256": tree_hash(entries), "n_bundles": len(entries),
            "decay_model": manifest["decay_model"], "exhad_commit": manifest["exhad_commit"],
            "exhad_variation": manifest["exhad_variation"],
            "pythia_version": manifest["pythia_version"],
            "n_templates": int(manifest["n_templates"]), "base_seed": int(manifest["base_seed"])}


def config_keys(variation) -> tuple:
    if variation.get("decay_templates"):
        return CONFIG_KEYS + ("decay_templates",)
    return CONFIG_KEYS


def expected_events_for_mass(m_S, n_pool) -> int:
    """Scalar events generated at ``m_S``: ``n_pool`` per OPEN b-hadron species.
    A species closes at ``m_S >= m_parent - m_recoil`` (Lambda_b at 4.50 GeV,
    below the ~4.78 GeV meson ceilings), so above 4.50 GeV only the three
    mesons contribute."""
    open_species = sum(1 for (_tag, m_parent, m_recoil, _frag) in production.B_SPECIES.values()
                       if m_S < m_parent - m_recoil)
    return n_pool * open_species


def campaign_config(variation, masses, n_pool, n_samples, seed) -> dict:
    payload = {
        "schema_version": 1,
        "variation": variation["name"],
        "axis": variation["axis"],
        "width_scheme": variation["width_scheme"],
        "grid_file": Path(variation["path"]).name,
        "grid_sha256": variation["sha256"],
        "n_parent_pool": int(n_pool),
        "n_scalar_events_expected": int(n_pool * len(production.B_SPECIES)),
        "n_decay_samples_per_hit": int(n_samples),
        "base_seed": int(seed),
        "parent_pool_seed": stable_seed(seed, variation["name"], "parent_pool"),
        "mass_grid_GeV": [float(m) for m in masses],
        "coupling_grid": {"log10_min": LOG_S2T_MIN, "log10_max": LOG_S2T_MAX, "n_points": N_S2T},
        "ray_backend": ray_backend_provenance(),
        "code_sha256": code_hashes(CODE_INPUTS),
    }
    if variation.get("decay_templates"):
        payload["decay_templates"] = templates_provenance(str(Path(variation["decay_templates"]).resolve()))
    payload["config_sha256"] = config_sha256(payload)
    return payload


def load_recorded_config(run_dir: Path, variation, masses, n_pool, n_samples, seed) -> dict:
    """A completed run's producer config, checked against today's campaign
    definition and against the producer commit it recorded (code identity is
    verified against that commit, not the collector's checkout)."""
    metadata = read_json(run_dir / "run_metadata.json")
    if metadata is None:
        raise RuntimeError(f"{run_dir.name}: missing retained run metadata")
    keys = config_keys(variation)
    missing = [key for key in keys if key not in metadata]
    if missing:
        raise RuntimeError(f"{run_dir.name}: incomplete retained config: {', '.join(missing)}")
    recorded = {key: metadata[key] for key in keys}
    if metadata.get("config_sha256") != config_sha256(recorded):
        raise RuntimeError(f"{run_dir.name}: retained config checksum mismatch")
    expected = campaign_config(variation, masses, n_pool, n_samples, seed)
    for key in keys:
        if key != "code_sha256" and recorded[key] != expected[key]:
            raise RuntimeError(f"{run_dir.name}: retained campaign definition differs at {key}")
    producer_commit = metadata.get("producer_git_head_at_start")
    if not producer_commit:
        raise RuntimeError(f"{run_dir.name}: missing producer commit")
    resolved = validate_recorded_code_state(producer_commit, recorded["code_sha256"], CODE_INPUTS)
    return {**recorded, "config_sha256": metadata["config_sha256"],
            "producer_git_head_at_start": resolved}


# --------------------------------------------------------------- run side --

@contextmanager
def bottom_grid(path: Path):
    old = os.environ.get(GRID_ENV)
    os.environ[GRID_ENV] = str(path)
    try:
        yield
    finally:
        if old is None:
            os.environ.pop(GRID_ENV, None)
        else:
            os.environ[GRID_ENV] = old


def ensure_run_metadata(run_dir: Path, config) -> None:
    path = run_dir / "run_metadata.json"
    old = read_json(path)
    if old is not None:
        if old.get("config_sha256") != config["config_sha256"]:
            raise RuntimeError(f"campaign configuration changed for {run_dir.name}; use a new "
                               "--scratch-dir or restore the recorded inputs")
        return
    atomic_json(path, {
        **config,
        "created_utc": utc_now(),
        "producer_git_head_at_start": git_head(),
        "independence_policy": ("fresh variation-specific FONLL parent pool, scalar four-vectors, "
                                "geometry, decay/reconstruction MC, and sensitivity scan; no reuse "
                                "or reweighting across variations"),
    })


def vector_paths(run_dir: Path, mass: float):
    vector = run_dir / "llp_4vectors" / f"mS_{production._mass_label(mass)}.csv"
    return vector, vector.with_suffix(".meta.json")


def valid_vector(vector: Path, metadata: Path, config, mass: float):
    meta = read_json(metadata) if vector.exists() else None
    if meta is None:
        return None
    if (meta.get("config_sha256") != config["config_sha256"]
            or float(meta.get("mass_GeV", -1)) != float(mass)
            or meta.get("bytes") != vector.stat().st_size
            or meta.get("sha256") != sha256_file(vector)):
        return None
    return meta


def write_vector(vector: Path, generated, config, variation, mass, seed):
    weight, energy, px, py, pz = generated
    matrix = np.column_stack([weight, energy, px, py, pz])
    vector.parent.mkdir(parents=True, exist_ok=True)
    tmp = vector.with_suffix(vector.suffix + f".tmp.{os.getpid()}")
    with open(tmp, "w") as fh:
        np.savetxt(fh, matrix, delimiter=",", fmt="%.8e")
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, vector)
    meta = {"config_sha256": config["config_sha256"], "variation": variation["name"],
            "grid_sha256": variation["sha256"], "mass_GeV": float(mass),
            "production_seed": int(seed), "rows": int(len(matrix)),
            "bytes": vector.stat().st_size, "sha256": sha256_file(vector)}
    atomic_json(vector.with_suffix(".meta.json"), meta)
    return meta


def result_path(run_dir: Path, mass: float) -> Path:
    return run_dir / "results" / f"mS_{production._mass_label(mass)}.json"


def valid_result(path: Path, config):
    result = read_json(path)
    return result if result is not None and result.get("config_sha256") == config["config_sha256"] else None


def raw_artifact_inventory(run_dir: Path, config, masses) -> dict:
    vectors, geometry, results = [], [], []
    for mass in masses:
        vector, meta_path = vector_paths(run_dir, mass)
        meta = valid_vector(vector, meta_path, config, mass)
        if meta is None:
            raise RuntimeError(f"{run_dir.name} m={mass:.3f}: vector failed final checksum validation")
        vectors.extend([file_record(vector, run_dir, known_sha=meta["sha256"]),
                        file_record(meta_path, run_dir)])
        geometry_path = run_dir / "geometry_cache" / f"geom_{vector.stem}.npz"
        if not geometry_path.exists():
            raise RuntimeError(f"{run_dir.name} m={mass:.3f}: missing geometry cache")
        geometry.append(file_record(geometry_path, run_dir))
        rp = result_path(run_dir, mass)
        if valid_result(rp, config) is None:
            raise RuntimeError(f"{run_dir.name} m={mass:.3f}: result failed final validation")
        results.append(file_record(rp, run_dir))
    return {"vectors": stage_record(vectors), "geometry": stage_record(geometry),
            "results": stage_record(results)}


def validate_retained_completion(run_dir: Path, config, masses, required_state="compacted"):
    completion = read_json(run_dir / "complete.json")
    if completion is None or completion.get("state") != required_state:
        return None
    if completion.get("config_sha256") != config["config_sha256"]:
        raise RuntimeError(f"{run_dir.name}: completion config checksum mismatch")
    if completion.get("mass_grid_GeV", []) != [float(m) for m in masses]:
        raise RuntimeError(f"{run_dir.name}: completion mass grid mismatch")
    curve_path = run_dir / completion["curve_file"]
    artifacts_path = run_dir / completion["artifacts_manifest_file"]
    if not curve_path.exists() or sha256_file(curve_path) != completion.get("curve_sha256"):
        raise RuntimeError(f"{run_dir.name}: retained curve checksum mismatch")
    if not artifacts_path.exists() or sha256_file(artifacts_path) != completion.get("artifacts_manifest_sha256"):
        raise RuntimeError(f"{run_dir.name}: artifact manifest checksum mismatch")
    artifacts = json.loads(artifacts_path.read_text())
    if artifacts.get("config_sha256") != config["config_sha256"]:
        raise RuntimeError(f"{run_dir.name}: artifact manifest config mismatch")
    entries = []
    for mass in masses:
        rp = result_path(run_dir, mass)
        if valid_result(rp, config) is None:
            raise RuntimeError(f"{run_dir.name} m={mass:.3f}: retained result is invalid; "
                               "use a new scratch directory to regenerate the compacted run")
        entries.append(file_record(rp, run_dir))
    if stage_record(entries)["tree_sha256"] != artifacts["stages"]["results"]["tree_sha256"]:
        raise RuntimeError(f"{run_dir.name}: retained result tree checksum mismatch")
    return completion, artifacts


def finish_compaction(run_dir: Path, completion) -> None:
    import shutil
    for dirname in completion["raw_directories"]:
        raw_dir = run_dir / dirname
        if raw_dir.exists():
            shutil.rmtree(raw_dir)
        if raw_dir.exists():
            raise RuntimeError(f"{run_dir.name}: failed to reclaim {raw_dir}")
    atomic_json(run_dir / "complete.json", {**completion, "state": "compacted", "compacted_utc": utc_now()})
    append_log(run_dir, "state=compacted raw vectors and geometry reclaimed")


def finalize_and_compact(run_dir: Path, variation, config, masses) -> None:
    rows = []
    for mass in masses:
        result = valid_result(result_path(run_dir, mass), config)
        if result is None:
            raise RuntimeError(f"{variation['name']} m={mass:.3f}: cannot finalize invalid result")
        rows.append(result)
    curve_path = run_dir / "curve.csv"
    atomic_csv(pd.DataFrame(rows).sort_values("mass_GeV"), curve_path)
    stages = raw_artifact_inventory(run_dir, config, masses)
    artifacts = {
        "schema_version": 1, "variation": variation["name"],
        "config_sha256": config["config_sha256"], "created_utc": utc_now(), "stages": stages,
        "regeneration": {"policy": "deterministic regeneration from pinned code, grid, and seeds",
                         "base_seed": config["base_seed"], "parent_pool_seed": config["parent_pool_seed"],
                         "code_sha256": config["code_sha256"], "grid_sha256": config["grid_sha256"],
                         "ray_backend": config["ray_backend"]},
    }
    artifacts_path = run_dir / "artifacts_manifest.json"
    atomic_json(artifacts_path, artifacts)
    completion = {
        "schema_version": 1, "state": "complete_raw", "variation": variation["name"],
        "config_sha256": config["config_sha256"], "completed_raw_utc": utc_now(),
        "n_masses": len(rows), "mass_grid_GeV": [float(m) for m in masses],
        "curve_file": curve_path.name, "curve_sha256": sha256_file(curve_path),
        "artifacts_manifest_file": artifacts_path.name,
        "artifacts_manifest_sha256": sha256_file(artifacts_path),
        "raw_directories": ["llp_4vectors", "geometry_cache"],
    }
    atomic_json(run_dir / "complete.json", completion)
    append_log(run_dir, f"state=complete_raw vectors_tree={stages['vectors']['tree_sha256']} "
                        f"geometry_tree={stages['geometry']['tree_sha256']} "
                        f"results_tree={stages['results']['tree_sha256']}")
    finish_compaction(run_dir, completion)


def run_variation(variation, scratch_dir: Path, masses, n_pool, n_samples, seed,
                  thresholds=(3.0,)):
    """Run one fully independent physics variation, resuming per mass."""
    from ...geometry.raycast import get_mesh
    scratch_dir = Path(scratch_dir)
    run_dir = scratch_dir / "runs" / variation["name"]
    config = campaign_config(variation, production.MASS_GRID, n_pool, n_samples, seed)
    with variation_lock(run_dir):
        ensure_run_metadata(run_dir, config)
        marker = read_json(run_dir / "complete.json")
        if marker is not None and marker.get("state") == "complete_raw":
            if validate_retained_completion(run_dir, config, masses, required_state="complete_raw") is None:
                raise RuntimeError(f"{variation['name']}: invalid complete_raw marker")
            finish_compaction(run_dir, marker)
            return variation["name"]
        if marker is not None and marker.get("state") == "compacted":
            validate_retained_completion(run_dir, config, masses)
            print(f"[{variation['name']}] already compacted", flush=True)
            return variation["name"]
        complete = sum(valid_result(result_path(run_dir, m), config) is not None for m in masses)
        print(f"[{variation['name']}] resume {complete}/{len(masses)} masses", flush=True)
        append_log(run_dir, f"resume results={complete}/{len(masses)}")
        if complete == len(masses):
            finalize_and_compact(run_dir, variation, config, masses)
            return variation["name"]

        pool_seed = config["parent_pool_seed"]
        with bottom_grid(variation["path"]):
            sigma = get_sigma_total("bottom")
            pool = sample_meson_4vectors(n_pool, "bottom", rng=np.random.default_rng(pool_seed))
        expected_sigma = variation["trapezoid_integral_pb"]
        # The manifest records a rounded decimal integral, while the parser
        # integrates the stored grid values at full float precision.
        if not np.isclose(sigma, expected_sigma, rtol=1e-8, atol=0.0):
            raise RuntimeError(f"{variation['name']}: grid integral {sigma} != manifest {expected_sigma}")

        paths = ModelPaths("bc4", vectors=run_dir / "llp_4vectors", templates=run_dir / "decay_templates",
                           analysis=run_dir, geometry=run_dir / "geometry_cache")
        spec = ScalarSpec(paths, templates_dir=variation.get("decay_templates"),
                          width_scheme=variation["width_scheme"])
        cfg = ScanConfig(decay_samples=n_samples, thresholds=tuple(thresholds))
        mesh = get_mesh()
        for index, mass in enumerate(masses):
            mass = float(mass)
            rp = result_path(run_dir, mass)
            if valid_result(rp, config) is not None:
                continue
            label = production._mass_label(mass)
            production_seed = stable_seed(seed, variation["name"], label, "production")
            reco_seed = stable_seed(seed, variation["name"], label, "reconstruction")
            vector, meta_path = vector_paths(run_dir, mass)
            meta = valid_vector(vector, meta_path, config, mass)
            if meta is None:
                generated = production.generate_scalar_4vectors(
                    mass, n_pool, np.random.default_rng(production_seed), sigma_bottom=sigma, pool=pool)
                meta = write_vector(vector, generated, config, variation, mass, production_seed)
            expected_events = expected_events_for_mass(mass, n_pool)
            if meta["rows"] != expected_events:
                raise RuntimeError(f"{variation['name']} m={mass}: expected {expected_events} scalar "
                                   f"events ({expected_events // n_pool} open species x {n_pool}), "
                                   f"got {meta['rows']}")
            result = run_point(spec, MassPoint(mass), cfg, mesh, rng=np.random.default_rng(reco_seed))
            if result is None:
                raise RuntimeError(f"{variation['name']} m={mass}: no production output "
                                   "(closed production, no hits, or no decay template)")
            record = {
                **jsonable(result.row),
                "variation": variation["name"], "axis": variation["axis"],
                "width_scheme": variation["width_scheme"], "mass_GeV": mass,
                "config_sha256": config["config_sha256"], "grid_sha256": variation["sha256"],
                "vector_sha256": meta["sha256"], "parent_pool_seed": int(pool_seed),
                "production_seed": int(production_seed), "reconstruction_seed": int(reco_seed),
                "n_parent_pool": int(n_pool), "n_decay_samples_per_hit": int(n_samples),
            }
            atomic_json(rp, record)
            append_log(run_dir, f"mass={mass:.3f} result={rp.name} vector_sha256={meta['sha256']} "
                                f"production_seed={production_seed} reconstruction_seed={reco_seed}")
            print(f"[{variation['name']}] {index + 1}/{len(masses)} m={mass:.3f} "
                  f"hits={record['n_hits']} peak={record['peak_N']:.3g} "
                  f"island=[{record['u2_min']}, {record['u2_max']}]", flush=True)
        finalize_and_compact(run_dir, variation, config, masses)
        return variation["name"]


# ----------------------------------------------------------- collect side --

def load_campaign(variations, scratch_dir, masses, n_pool, n_samples, seed):
    """The raw variation table and per-run provenance of a compacted campaign."""
    scratch_dir = Path(scratch_dir)
    frames, run_records, missing = [], [], []
    for variation in variations:
        run_dir = scratch_dir / "runs" / variation["name"]
        config = load_recorded_config(run_dir, variation, production.MASS_GRID, n_pool, n_samples, seed)
        completed = validate_retained_completion(run_dir, config, masses)
        if completed is None:
            missing.append(f"{variation['name']}:not_compacted")
            continue
        completion, artifacts = completed
        rows = []
        for mass in masses:
            result = valid_result(result_path(run_dir, mass), config)
            if result is None:
                missing.append(f"{variation['name']}:{mass:.3f}")
            else:
                rows.append(result)
        if len(rows) == len(masses):
            frames.append(pd.DataFrame(rows).sort_values("mass_GeV"))
            run_records.append({
                "variation": variation["name"], "config_sha256": config["config_sha256"],
                "n_masses": len(rows), "parent_pool_seed": config["parent_pool_seed"],
                "producer_git_head_at_start": config["producer_git_head_at_start"],
                "producer_code_sha256": config["code_sha256"],
                "completion_state": completion["state"], "curve_sha256": completion["curve_sha256"],
                "artifacts_manifest_sha256": completion["artifacts_manifest_sha256"],
                "vector_tree_sha256": artifacts["stages"]["vectors"]["tree_sha256"],
                "geometry_tree_sha256": artifacts["stages"]["geometry"]["tree_sha256"],
                "results_tree_sha256": artifacts["stages"]["results"]["tree_sha256"],
            })
    if missing:
        raise RuntimeError(f"campaign incomplete: {len(missing)} missing mass results; "
                           f"first: {', '.join(missing[:20])}")
    raw = pd.concat(frames, ignore_index=True).sort_values(["axis", "variation", "mass_GeV"])
    return raw, run_records


def _control_summary(band, boundary):
    prefix = f"{boundary}_numerical_repeat"
    flagged = band[f"{prefix}_not_subdominant"].eq(True)
    summary = {"n_masses_not_subdominant": int(flagged.sum()),
               "masses_not_subdominant_GeV": [float(v) for v in band.loc[flagged, "mass_GeV"]],
               "max_abs_dex": float(band[f"{prefix}_max_abs_dex"].max()),
               "max_fractional": float(band[f"{prefix}_max_fractional"].max())}
    for component in ("scale", "pdf", "bottom_mass", "decay_model"):
        values = band[f"{prefix}_to_{component}_ratio"].replace([np.inf, -np.inf], np.nan)
        summary[f"max_ratio_to_{component}"] = float(values.max())
    return summary


def _topology_summary(band, prefix):
    out = {}
    for reference in ("campaign", "canonical"):
        flag = band[f"{prefix}_topology_differs_from_{reference}"].fillna(False).astype(bool)
        column = f"{prefix}_topology_difference_variations_from_{reference}"
        out[f"topology_vs_{reference}"] = {
            "n_differences": int(flag.sum()),
            "difference_masses_GeV": [float(v) for v in band.loc[flag, "mass_GeV"]],
            "difference_variations_by_mass": {
                f"{float(row.mass_GeV):g}": getattr(row, column)
                for row in band.loc[flag, ["mass_GeV", column]].itertuples(index=False)},
        }
    return out


def collect_campaign(variations, scratch_dir, masses, n_pool, n_samples, seed,
                     central_curve, curves_out, band_out, manifest_out):
    """Require compacted runs and publish the raw table, the envelope and its provenance."""
    raw, run_records = load_campaign(variations, scratch_dir, masses, n_pool, n_samples, seed)
    band = combine_band(raw, central_curve)
    atomic_csv(raw, Path(curves_out))
    atomic_csv(band, Path(band_out))
    source_manifest = Path(variations[0]["path"]).parent / "variation_manifest.json"
    producer_states = {json.dumps({"git_head": r["producer_git_head_at_start"],
                                   "code_sha256": r["producer_code_sha256"]}, sort_keys=True)
                       for r in run_records}
    if len(producer_states) != 1:
        raise RuntimeError("campaign variations were produced from different code states")
    producer_state = json.loads(next(iter(producer_states)))
    manifest = {
        "artifact": "GRENDEL BC4 independent full-statistics uncertainty campaign",
        "schema_version": 2,
        "published_utc": utc_now(),
        "producer": {"git_head": producer_state["git_head"], "code_sha256": producer_state["code_sha256"]},
        "collector": {"git_head": git_head()},
        "independence_policy": ("Every variation has a fresh FONLL parent sample, scalar four-vectors, "
                                "geometry, decay/reconstruction MC, and sensitivity scan. No event, "
                                "weight, geometry, or detector-outcome reuse/reweighting across variations."),
        "campaign": {
            "n_fonll_variations": sum(v["axis"] in ("central", "scale", "pdf", "mass") for v in variations),
            "n_decay_model_variations": sum(v["axis"] == "decay_model" for v in variations),
            "n_numerical_control_variations": len(NUMERICAL_CONTROL_VARIATIONS),
            "n_total_variations": len(variations), "n_masses": len(masses),
            "n_parent_pool_per_variation": n_pool,
            "n_scalar_events_per_mass": n_pool * len(production.B_SPECIES),
            "n_decay_samples_per_hit": n_samples, "base_seed": seed,
            "mass_grid_GeV": [float(m) for m in masses],
            "source_fonll_manifest_sha256": sha256_file(source_manifest),
        },
        "combination": {
            "space": "log10(sin^2 theta)",
            "scale": "asymmetric envelope of central plus six coherent scale grids",
            "pdf": "16th and 84th percentiles (linear quantiles) of 100 NNPDF replica boundaries; "
                   "sample standard deviation retained for audit",
            "bottom_mass": "extrema of central plus mb=4.5/5.0 GeV variations",
            "decay_model": "interval between central and the independent central-FONLL decay-model "
                           "alternates (LO-ChPT widths below 2 GeV and perturbative spectator widths "
                           "above; exHad when supplemented)",
            "display": "single_source_variation_envelope: outermost endpoint among the scale, "
                       "PDF-percentile, bottom-mass, and decay-model intervals",
            "interpretation": "one source varied at a time; no quadrature, no simultaneous-source "
                              "coverage, and not a confidence band",
            "numerical_controls": "two fresh-seed same-physics central repeats, reported separately "
                                  "and excluded from all physics intervals and the display envelope",
            "rebase": "component dex shifts applied to the canonical central curve",
        },
        "numerical_control_summary": {
            "definition": "absolute boundary shifts of two fresh-seed, same-physics central repeats "
                          "relative to the independent campaign central",
            "envelope_policy": "excluded from every physics interval and display envelope",
            "subdominant_criterion": "repeat maximum absolute dex shift is smaller than the largest "
                                     "one-source physical interval displacement at that mass/boundary",
            "u2_min": _control_summary(band, "u2_min"), "u2_max": _control_summary(band, "u2_max"),
            **_topology_summary(band, "numerical_control"),
        },
        "physical_variation_topology_summary": _topology_summary(band, "physical_variation"),
        "variations": [{"name": v["name"], "axis": v["axis"], "width_scheme": v["width_scheme"],
                        "grid_file": Path(v["path"]).name, "grid_sha256": v["sha256"],
                        "trapezoid_integral_pb": v["trapezoid_integral_pb"],
                        **next(r for r in run_records if r["variation"] == v["name"])}
                       for v in variations],
        "outputs": {Path(curves_out).name: {"sha256": sha256_file(Path(curves_out)), "rows": len(raw)},
                    Path(band_out).name: {"sha256": sha256_file(Path(band_out)), "rows": len(band)},
                    Path(central_curve).name: {"sha256": sha256_file(Path(central_curve)),
                                               "role": "canonical central"}},
        "software": {**software(), "ray_backend": ray_backend_provenance()},
        "limitations": ["No FONLL alpha_s companion grids are included.",
                        "The decay alternate is a model envelope, not a Gaussian error.",
                        "The display envelope is not a confidence interval or simultaneous-source band.",
                        "Detector-response and background systematics are outside this theory envelope."],
    }
    atomic_json(Path(manifest_out), manifest)
    return raw, band, manifest


def campaign_status(variations, scratch_dir, masses, n_pool, n_samples, seed):
    total = len(variations) * len(masses)
    done = 0
    for variation in variations:
        config = campaign_config(variation, production.MASS_GRID, n_pool, n_samples, seed)
        run_dir = Path(scratch_dir) / "runs" / variation["name"]
        count = sum(valid_result(result_path(run_dir, m), config) is not None for m in masses)
        marker = read_json(run_dir / "complete.json")
        state = marker.get("state", "in_progress") if marker else "in_progress"
        done += count
        print(f"{variation['name']:28s} {count:3d}/{len(masses)} {state}")
    print(f"TOTAL {done}/{total} mass-variation results ({100 * done / max(total, 1):.2f}%)")
    return done, total


# --------------------------------------------------------------------- CLI --

def add_common_options(parser):
    parser.add_argument("--grid-dir", type=Path, default=os.environ.get("GRENDEL_FONLL_GRID_DIR"),
                        required=os.environ.get("GRENDEL_FONLL_GRID_DIR") is None,
                        help="the FONLL variation campaign (variation_manifest.json + grids)")
    parser.add_argument("--scratch-dir", type=Path, default=os.environ.get("GRENDEL_BC4_CAMPAIGN_DIR"),
                        required=os.environ.get("GRENDEL_BC4_CAMPAIGN_DIR") is None,
                        help="campaign workspace for runs/<variation>/")
    parser.add_argument("--n-pool", type=int, default=production.N_POOL_DEFAULT,
                        help="bottom parent events; the four b-hadron species give 4*n-pool scalars per mass")
    parser.add_argument("--n-samples", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--masses", type=float, nargs="+", default=None)
    parser.add_argument("--variations", nargs="+", default=None,
                        help="variation names or axes; default complete campaign")
    parser.add_argument("--thresholds", nargs="+", type=float, default=[3.0])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    run_parser = sub.add_parser("run", help="run/resume selected independent variations")
    add_common_options(run_parser)
    run_parser.add_argument("--workers", type=int, default=1, help="independent variation processes")
    status_parser = sub.add_parser("status", help="report completed mass checkpoints")
    add_common_options(status_parser)
    collect_parser = sub.add_parser("collect", help="require completeness and publish the band")
    add_common_options(collect_parser)
    collect_parser.add_argument("--central-curve", type=Path, required=True,
                                help="the canonical BC4 sensitivity.csv the envelope is drawn around")
    collect_parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    variations = select_variations(discover_variations(args.grid_dir.resolve()), args.variations)
    masses = [float(m) for m in (args.masses or production.MASS_GRID)]
    if args.command == "status":
        campaign_status(variations, args.scratch_dir, masses, args.n_pool, args.n_samples, args.seed)
        return 0
    if args.command == "collect":
        if args.masses or args.variations:
            raise ValueError("collection requires the complete mass/variation campaign")
        out = args.out_dir
        out.mkdir(parents=True, exist_ok=True)
        raw, band, _ = collect_campaign(
            variations, args.scratch_dir, masses, args.n_pool, args.n_samples, args.seed,
            args.central_curve, out / "bc4_uncertainty_variations.csv",
            out / "bc4_single_source_variation_envelope.csv", out / "UNCERTAINTY_MANIFEST.json")
        print(f"wrote {len(raw)} variation rows and {len(band)} band rows -> {out}")
        return 0
    if args.workers < 1:
        raise ValueError("--workers must be positive")
    args.scratch_dir.mkdir(parents=True, exist_ok=True)
    thresholds = tuple(args.thresholds)
    if args.workers == 1:
        for variation in variations:
            run_variation(variation, args.scratch_dir, masses, args.n_pool, args.n_samples, args.seed, thresholds)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(run_variation, v, args.scratch_dir, masses, args.n_pool,
                                       args.n_samples, args.seed, thresholds): v["name"]
                       for v in variations}
            for future in as_completed(futures):
                name = futures[future]
                try:
                    future.result()
                    print(f"[orchestrator] complete: {name}", flush=True)
                except Exception as exc:
                    print(f"[orchestrator] FAILED: {name}: {exc}", file=sys.stderr, flush=True)
                    raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
