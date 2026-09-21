"""Add the exHad decay model to the BC4 variation envelope as one more independent variation.

A campaign spans its decay-model axis with two points, the Winkler central
and the LO-ChPT/spectator alternate. exHad's ``scalar-1809`` widths lie
outside that interval above 2 GeV, so the axis needs a third point. This
module runs that point under the campaign's independence policy -- its own
FONLL parent pool (central grid), scalar four-vectors, ray casts,
decay/reconstruction Monte Carlo and coupling scan at the campaign's pool
size and seed policy -- through the very same ``campaign.run_variation``,
and combines it with the campaign's raw curves into a widened envelope
written next to them:

    <out>/bc4_uncertainty_variations_exhad.csv            the new variation's rows
    <out>/bc4_single_source_variation_envelope_exhad.csv  the widened envelope
    <out>/UNCERTAINTY_MANIFEST_exhad.json                 provenance

The decay templates are a fresh exHad sample too (their own seed, derived
from the variation name), so the decay stage shares nothing with the exHad
central-curve run.

    python -m grendel.band.scalar.exhad_variation templates --exhad-python .../bin/python
    python -m grendel.band.scalar.exhad_variation run
    python -m grendel.band.scalar.exhad_variation status
    python -m grendel.band.scalar.exhad_variation collect --campaign-raw ... --campaign-band ... \\
        --campaign-manifest ... --central-curve ... --out-dir ...
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from ...io.atomic import atomic_csv, atomic_json, sha256_file
from ...io.paths import repo_root
from ...models.scalar import production
from ..campaign_store import code_hashes, git_head, software, stable_seed, utc_now
from . import campaign as camp
from .combine import combine_band
from .variations import EXHAD_DECAY_VARIATION, exhad_variation

TEMPLATES_DIRNAME = "templates_exhad_1809"
EXHAD_MODEL = "scalar-1809"


def campaign_parameters(manifest_path: Path) -> dict:
    """Pool size, decay samples, seed and mass grid of the campaign, read from
    its manifest so the new variation cannot silently run at another size."""
    manifest = json.loads(Path(manifest_path).read_text())
    c = manifest["campaign"]
    grid = [float(m) for m in c["mass_grid_GeV"]]
    if grid != [float(m) for m in production.MASS_GRID]:
        raise RuntimeError("the BC4 mass grid differs from the campaign grid")
    return {"n_pool": int(c["n_parent_pool_per_variation"]), "n_samples": int(c["n_decay_samples_per_hit"]),
            "seed": int(c["base_seed"]), "masses": grid, "manifest_sha256": sha256_file(manifest_path)}


def templates_seed(seed: int) -> int:
    """Base seed of the variation's own exHad template sample."""
    return stable_seed(seed, EXHAD_DECAY_VARIATION, "templates")


def templates_dir(scratch: Path) -> Path:
    return Path(scratch) / TEMPLATES_DIRNAME


def cmd_templates(args) -> int:
    out = templates_dir(args.scratch_dir)
    params = campaign_parameters(args.campaign_manifest)
    seed = templates_seed(params["seed"])
    if (out / "MANIFEST.json").exists() and not args.force:
        print(f"templates already generated under {out} (use --force to regenerate)")
        return 0
    command = [str(args.exhad_python), "-m", "grendel.models.scalar.templates_exhad", "--model", EXHAD_MODEL,
               "--n-templates", str(args.n_templates), "--seed", str(seed),
               "--workers", str(args.workers), "--out", str(out)]
    print("$", " ".join(command), flush=True)
    return subprocess.run(command, cwd=repo_root(), env={**os.environ, "PYTHONUNBUFFERED": "1"}).returncode


def _check_templates(variation, masses) -> None:
    root = Path(variation["decay_templates"])
    missing = [m for m in masses if not (root / f"templates_{production._mass_label(m)}.npz").exists()]
    if missing:
        raise FileNotFoundError(f"{len(missing)} template bundles missing under {root}; first {missing[:3]}")
    manifest = json.loads((root / "MANIFEST.json").read_text())
    if manifest["decay_model"] != f"exhad:{EXHAD_MODEL}":
        raise RuntimeError(f"unexpected template model {manifest['decay_model']!r}")


def cmd_run(args) -> int:
    params = campaign_parameters(args.campaign_manifest)
    variation = exhad_variation(args.grid_dir, templates_dir(args.scratch_dir))
    _check_templates(variation, params["masses"])
    manifest = json.loads((Path(variation["decay_templates"]) / "MANIFEST.json").read_text())
    if int(manifest["base_seed"]) != templates_seed(params["seed"]):
        raise RuntimeError(f"template seed {manifest['base_seed']} is not the variation's "
                           f"{templates_seed(params['seed'])}")
    Path(args.scratch_dir).mkdir(parents=True, exist_ok=True)
    print(f"[{variation['name']}] n_pool={params['n_pool']} n_samples={params['n_samples']} "
          f"seed={params['seed']} masses={len(params['masses'])} templates={variation['decay_templates']}",
          flush=True)
    camp.run_variation(variation, args.scratch_dir, params["masses"], params["n_pool"],
                       params["n_samples"], params["seed"])
    return 0


def cmd_status(args) -> int:
    params = campaign_parameters(args.campaign_manifest)
    variation = exhad_variation(args.grid_dir, templates_dir(args.scratch_dir))
    camp.campaign_status([variation], args.scratch_dir, params["masses"], params["n_pool"],
                         params["n_samples"], params["seed"])
    return 0


def containment(curve: pd.DataFrame, band: pd.DataFrame) -> dict:
    """How much of a curve lies inside an envelope, per edge (dex breaches)."""
    curve = curve.rename(columns={"has_sensitivity": "has_sensitivity_curve",
                                  "u2_min": "u2_min_curve", "u2_max": "u2_max_curve"})
    merged = band.merge(curve, on="mass_GeV")
    merged = merged[merged["has_sensitivity"].astype(bool) & merged["has_sensitivity_curve"].astype(bool)]
    out = {}
    for edge in ("u2_min", "u2_max"):
        value = merged[f"{edge}_curve"].to_numpy(float)
        lo = merged[f"{edge}_envelope_lo"].to_numpy(float)
        hi = merged[f"{edge}_envelope_hi"].to_numpy(float)
        ok = np.isfinite(value) & np.isfinite(lo) & np.isfinite(hi)
        value, lo, hi = value[ok], lo[ok], hi[ok]
        inside = (value >= lo) & (value <= hi)
        excess = np.where(value > hi, np.log10(value / hi), np.log10(value / lo))
        excess = np.where(inside, 0.0, excess)
        out[edge] = {"n": int(len(inside)), "inside": int(inside.sum()),
                     "worst_breach_dex": float(np.abs(excess).max()) if len(excess) else None}
    return out


def reproduce_campaign_envelope(raw: pd.DataFrame, central_curve: Path, band_path: Path) -> dict:
    """Guard: the collector must still reproduce the campaign envelope from
    the campaign's raw curves before it is trusted to extend it."""
    recomputed = combine_band(raw, central_curve)
    published = pd.read_csv(band_path)
    worst = 0.0
    for column in published.columns:
        if column not in recomputed:
            raise RuntimeError(f"collector no longer produces column {column}")
        x, y = recomputed[column], published[column]
        if pd.api.types.is_numeric_dtype(y) and pd.api.types.is_numeric_dtype(x):
            xv, yv = x.to_numpy(float), y.to_numpy(float)
            if not np.array_equal(np.isfinite(xv), np.isfinite(yv)):
                raise RuntimeError(f"column {column} not reproduced (NaN pattern)")
            finite = np.isfinite(xv)
            if finite.any():
                worst = max(worst, float(np.abs(xv[finite] - yv[finite]).max()))
                if not np.allclose(xv[finite], yv[finite], rtol=1e-9, atol=1e-9):
                    raise RuntimeError(f"column {column} not reproduced")
        elif not (x.fillna("").astype(str).to_numpy() == y.fillna("").astype(str).to_numpy()).all():
            raise RuntimeError(f"column {column} not reproduced")
    return {"reproduced": True, "n_columns": int(len(published.columns)), "max_abs_difference": worst}


def cmd_collect(args) -> int:
    params = campaign_parameters(args.campaign_manifest)
    variation = exhad_variation(args.grid_dir, templates_dir(args.scratch_dir))
    run_dir = Path(args.scratch_dir) / "runs" / variation["name"]
    config = camp.load_recorded_config(run_dir, variation, params["masses"], params["n_pool"],
                                       params["n_samples"], params["seed"])
    completed = camp.validate_retained_completion(run_dir, config, params["masses"])
    if completed is None:
        raise RuntimeError(f"{variation['name']} is not compacted; run `run` to completion first")
    completion, artifacts = completed
    rows = []
    for mass in params["masses"]:
        result = camp.valid_result(camp.result_path(run_dir, mass), config)
        if result is None:
            raise RuntimeError(f"{variation['name']} m={mass:.3f}: retained result invalid")
        rows.append(result)
    new_rows = pd.DataFrame(rows).sort_values("mass_GeV").reset_index(drop=True)
    campaign_raw = pd.read_csv(args.campaign_raw)
    campaign_manifest = json.loads(Path(args.campaign_manifest).read_text())
    recorded = campaign_manifest["outputs"][Path(args.campaign_raw).name]["sha256"]
    if sha256_file(args.campaign_raw) != recorded:
        raise RuntimeError("the campaign raw variation curves do not match their manifest hash")
    pools = set(campaign_raw["n_parent_pool"].astype(int).unique())
    if pools != {params["n_pool"]} or set(new_rows["n_parent_pool"].astype(int)) != {params["n_pool"]}:
        raise RuntimeError(f"pool sizes differ: campaign {pools}, new {set(new_rows['n_parent_pool'])}")
    if EXHAD_DECAY_VARIATION in set(campaign_raw["variation"]):
        raise RuntimeError("the campaign raw curves already contain the exHad variation")
    reproduction = reproduce_campaign_envelope(campaign_raw, args.central_curve, args.campaign_band)
    ordered = list(campaign_raw.columns) + [c for c in new_rows.columns if c not in campaign_raw.columns]
    new_rows = new_rows.reindex(columns=ordered)
    raw = pd.concat([campaign_raw, new_rows], ignore_index=True).sort_values(["axis", "variation", "mass_GeV"])
    band = combine_band(raw, args.central_curve)
    campaign_band = pd.read_csv(args.campaign_band)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_raw = out_dir / "bc4_uncertainty_variations_exhad.csv"
    out_band = out_dir / "bc4_single_source_variation_envelope_exhad.csv"
    out_manifest = out_dir / "UNCERTAINTY_MANIFEST_exhad.json"
    atomic_csv(new_rows, out_raw)
    atomic_csv(band, out_band)

    def widening(edge):
        lo_new, hi_new = band[f"{edge}_decay_model_envelope_lo"], band[f"{edge}_decay_model_envelope_hi"]
        lo_old, hi_old = campaign_band[f"{edge}_decay_model_envelope_lo"], campaign_band[f"{edge}_decay_model_envelope_hi"]
        down, up = np.log10(lo_old / lo_new), np.log10(hi_new / hi_old)
        fd, fu = np.isfinite(down), np.isfinite(up)
        shift = band[f"{edge}_decay_model_{EXHAD_DECAY_VARIATION}_shift_dex"]
        return {"n_masses_with_interval": int(fd.sum()),
                "n_masses_extended_down": int((down[fd] > 1e-12).sum()),
                "n_masses_extended_up": int((up[fu] > 1e-12).sum()),
                "max_extension_down_dex": float(down[fd].max()) if fd.any() else None,
                "max_extension_up_dex": float(up[fu].max()) if fu.any() else None,
                "exhad_shift_dex_vs_campaign_central": {
                    "n": int(shift.notna().sum()),
                    "median": float(shift.median()) if shift.notna().any() else None,
                    "min": float(shift.min()) if shift.notna().any() else None,
                    "max": float(shift.max()) if shift.notna().any() else None}}

    exhad_curve = pd.read_csv(args.exhad_central) if args.exhad_central and Path(args.exhad_central).exists() else None
    manifest = {
        "artifact": "GRENDEL BC4 uncertainty campaign, exHad decay-model supplement",
        "schema_version": 1,
        "published_utc": utc_now(),
        "status": "supplement: the campaign files are unchanged; this envelope adds the exHad scalar-1809 "
                  "decay model as a third point on the decay-model axis, run as one more independent variation",
        "independence_policy": campaign_manifest["independence_policy"],
        "extends": {Path(args.campaign_manifest).name: params["manifest_sha256"],
                    Path(args.campaign_raw).name: recorded,
                    Path(args.campaign_band).name: sha256_file(args.campaign_band),
                    Path(args.central_curve).name: sha256_file(args.central_curve),
                    "campaign_envelope_reproduced_by_this_collector": reproduction},
        "producer": {"git_head": git_head(), "code_sha256": code_hashes(camp.CODE_INPUTS),
                     "producer_git_head_at_start": config["producer_git_head_at_start"],
                     "code_sha256_at_start": config["code_sha256"]},
        "campaign": {"n_parent_pool": params["n_pool"],
                     "n_scalar_events_per_mass": params["n_pool"] * len(production.B_SPECIES),
                     "n_decay_samples_per_hit": params["n_samples"], "base_seed": params["seed"],
                     "n_masses": len(params["masses"]), "mass_grid_GeV": params["masses"]},
        "variation": {"name": variation["name"], "axis": variation["axis"], "width_scheme": variation["width_scheme"],
                      "grid_file": Path(variation["path"]).name, "grid_sha256": variation["sha256"],
                      "trapezoid_integral_pb": variation["trapezoid_integral_pb"],
                      "config_sha256": config["config_sha256"], "parent_pool_seed": config["parent_pool_seed"],
                      "decay_templates": config["decay_templates"], "completion_state": completion["state"],
                      "curve_sha256": completion["curve_sha256"],
                      "artifacts_manifest_sha256": completion["artifacts_manifest_sha256"],
                      "vector_tree_sha256": artifacts["stages"]["vectors"]["tree_sha256"],
                      "geometry_tree_sha256": artifacts["stages"]["geometry"]["tree_sha256"],
                      "results_tree_sha256": artifacts["stages"]["results"]["tree_sha256"],
                      "decay_model": "exHad (arXiv:2609.16104) scalar-1809 final states and lifetime: "
                                     "arXiv:1809.01876 widths with the corrections of arXiv:1904.10447, "
                                     "drawn from a fresh template sample per mass"},
        "combination": {**campaign_manifest["combination"],
                        "decay_model": "interval spanned by the central and the independent central-FONLL "
                                       "alternates: LO-ChPT widths below 2 GeV with perturbative spectator "
                                       "widths above (decay_chpt_spectator) and exHad scalar-1809 (decay_exhad_1809)"},
        "decay_model_axis_widening": {edge: widening(edge) for edge in ("u2_min", "u2_max")},
        "exhad_central_curve_containment": (
            {"curve": Path(args.exhad_central).name, "sha256": sha256_file(args.exhad_central),
             "campaign_envelope": containment(exhad_curve, campaign_band),
             "supplemented_envelope": containment(exhad_curve, band)} if exhad_curve is not None else None),
        "outputs": {out_raw.name: {"sha256": sha256_file(out_raw), "rows": int(len(new_rows))},
                    out_band.name: {"sha256": sha256_file(out_band), "rows": int(len(band))}},
        "software": {**software(), "ray_backend": camp.ray_backend_provenance()},
        "limitations": campaign_manifest["limitations"],
    }
    atomic_json(out_manifest, manifest)
    print(f"wrote {out_raw} ({len(new_rows)} rows), {out_band} ({len(band)} rows), {out_manifest}")
    print(json.dumps({"decay_model_axis_widening": manifest["decay_model_axis_widening"],
                      "exhad_central_curve_containment": manifest["exhad_central_curve_containment"]}, indent=1))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grid-dir", type=Path, default=os.environ.get("GRENDEL_FONLL_GRID_DIR"),
                    required=os.environ.get("GRENDEL_FONLL_GRID_DIR") is None)
    ap.add_argument("--scratch-dir", type=Path, default=os.environ.get("GRENDEL_BC4_EXHAD_CAMPAIGN_DIR"),
                    required=os.environ.get("GRENDEL_BC4_EXHAD_CAMPAIGN_DIR") is None,
                    help="workspace of this variation (its templates and run)")
    ap.add_argument("--campaign-manifest", type=Path, required=True,
                    help="UNCERTAINTY_MANIFEST.json of the campaign being extended")
    sub = ap.add_subparsers(dest="command", required=True)
    t = sub.add_parser("templates", help="generate the variation's own exHad template sample")
    t.add_argument("--exhad-python", type=Path, default=os.environ.get("GRENDEL_EXHAD_PYTHON"),
                   required=os.environ.get("GRENDEL_EXHAD_PYTHON") is None,
                   help="the interpreter of exHad's virtual environment")
    t.add_argument("--n-templates", type=int, default=20_000)
    t.add_argument("--workers", type=int, default=6)
    t.add_argument("--force", action="store_true")
    sub.add_parser("run", help="run/resume the exHad decay-model variation")
    sub.add_parser("status", help="report completed mass checkpoints")
    c = sub.add_parser("collect", help="combine with the campaign curves into the supplement")
    c.add_argument("--campaign-raw", type=Path, required=True, help="bc4_uncertainty_variations.csv")
    c.add_argument("--campaign-band", type=Path, required=True, help="bc4_single_source_variation_envelope.csv")
    c.add_argument("--central-curve", type=Path, required=True, help="the canonical BC4 sensitivity.csv")
    c.add_argument("--exhad-central", type=Path, default=None,
                   help="the exHad-decayed BC4 curve, for the containment check")
    c.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)
    return {"templates": cmd_templates, "run": cmd_run, "status": cmd_status, "collect": cmd_collect}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
