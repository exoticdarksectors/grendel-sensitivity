"""Drive the FONLL production-uncertainty band for the HNL exclusion curves.

Runs the production + scan chain once per *coherent* FONLL variation grid,
reusing the FONLL-independent channels, and leaves one ``sensitivity.csv``
per variation for ``grendel.band.quadrature`` to turn into boundary ribbons.

Coherence rules:

* **scale / PDF** -- the SAME variation grid is applied to BOTH bottom and
  charm. A scale choice or a proton-PDF replica is correlated across the two
  heavy flavours, so they move together.
* **mass** -- bottom and charm masses are INDEPENDENT, so each mass
  variation moves one quark's grid while the other stays central.

FONLL-dependent channels, regenerated per variation: ``bottom`` (Bmeson),
``charm`` (Dmeson), ``baryon`` (Bbaryon), ``induced_tau``. FONLL-independent
channels, reused from the central run via hardlinks: ``Bc``, ``Kmeson``,
``tau``, ``WZ`` (Bc borrows only the bottom shape with its own
normalisation, so it is held central).

Each variation runs in its own directory named after a hash of the two
grids' checksums and of the campaign configuration, so a directory can
never serve a geometry cache built from a different grid.

    python -m grendel.band.hnl.campaign --grid-dir G --campaign-dir C --flavor Ue Umu Utau
    python -m grendel.band.hnl.campaign --grid-dir G --campaign-dir C --variations central scale \\
        --flavor Umu --mass 1.0 2.0 --fonll-channels-only          # plumbing test
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

from ...io.paths import repo_root
from ...io.vectors import format_mass_for_filename
from ...models.hnl.mass_grid import MASS_GRID
from ...reco.acceptance import P_CUT
from ..campaign_store import git_head
from ..fonll_variations import grid_path, load_manifest, manifest_index

# run_all channel names regenerated per variation (FONLL-dependent)
FONLL_CHANNELS = ["bottom", "charm", "baryon", "induced_tau"]
# combine-channel label dirs reused from central (FONLL-independent)
REUSE_LABELS = ["Bc", "Kmeson", "tau", "WZ"]
# bottom/charm mass-variation grid tags
MASS_TAGS = {
    "mb_dn": ("mass_4p5", None),   # m_b = 4.50, charm central
    "mb_up": ("mass_5", None),     # m_b = 5.00, charm central
    "mc_dn": (None, "mass_1p3"),   # m_c = 1.30, bottom central
    "mc_up": (None, "mass_1p7"),   # m_c = 1.70, bottom central
}


def discover_variations(grid_dir: Path, which=None) -> list[dict]:
    """The coherent variation list from the campaign manifest."""
    by_key = manifest_index(load_manifest(grid_dir))

    def entry(tag, quark):
        e = by_key.get((quark, tag))
        if e is None:
            raise FileNotFoundError(f"manifest has no {quark} grid tagged {tag!r}")
        return e

    tags_present = {t for (_, t) in by_key}
    cb, cc = entry("central", "bottom"), entry("central", "charm")
    variations = [dict(name="central", axis="central", direction="0", bottom=cb, charm=cc)]
    for tag in sorted(t for t in tags_present if t.startswith("scale_")):
        variations.append(dict(name=tag, axis="scale", direction=tag,
                               bottom=entry(tag, "bottom"), charm=entry(tag, "charm")))
    for tag in sorted(t for t in tags_present if t.startswith("pdf_")):
        member = int(tag.split("_")[1])
        if member < 1:
            continue
        variations.append(dict(name=tag, axis="pdf", direction=str(member),
                               bottom=entry(tag, "bottom"), charm=entry(tag, "charm")))
    for name, (btag, ctag) in MASS_TAGS.items():
        if (btag and btag not in tags_present) or (ctag and ctag not in tags_present):
            continue
        variations.append(dict(name=name, axis="mass", direction=name,
                               bottom=entry(btag, "bottom") if btag else cb,
                               charm=entry(ctag, "charm") if ctag else cc))
    if which:
        sel = set(which)
        variations = [v for v in variations if v["name"] in sel or v["axis"] in sel]
    return variations


def run_tag(v: dict, campaign_config: dict | None = None) -> str:
    h = hashlib.sha256((v["bottom"]["sha256"] + v["charm"]["sha256"]).encode()).hexdigest()[:8]
    tag = f"band_{v['name']}_{h}"
    if campaign_config:
        payload = json.dumps(campaign_config, sort_keys=True, separators=(",", ":"))
        tag += f"_cfg{hashlib.sha256(payload.encode()).hexdigest()[:8]}"
    return tag


def code_revision() -> str:
    """The exact revision, marking uncommitted code explicitly."""
    root = repo_root()
    revision = git_head(root)
    dirty = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", "grendel"], cwd=root,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode != 0
    return f"{revision}-dirty" if dirty else revision


def campaign_configuration(args, flavors, masses) -> dict:
    return {"code_revision": code_revision(), "flavors": list(flavors),
            "masses_GeV": [float(m) for m in masses], "n_pool": int(args.n_pool),
            "production_seed": int(args.seed), "decay_samples": int(args.decay_samples),
            "max_hit_events": args.max_hit_events, "event_chunk": args.event_chunk,
            "analysis_seed_salt": args.analysis_seed_salt, "track_momentum_cut_GeV": float(P_CUT),
            "thresholds": [float(t) for t in args.thresholds], "reuse_from": args.reuse_from,
            "fonll_channels_only": bool(args.fonll_channels_only), "wz_nb_core": int(args.wz_nb_core)}


def run_dirs(campaign_dir: Path, tag: str):
    run_dir = campaign_dir / "runs" / tag
    return run_dir, run_dir / "llp_4vectors", run_dir / "analysis"


def environment(v: dict, grid_dir: Path, campaign_dir: Path, tag: str, templates_dir: Path) -> dict:
    run_dir, vec, analysis = run_dirs(campaign_dir, tag)
    env = os.environ.copy()
    env.update({
        "GRENDEL_WORK_DIR": str(run_dir),
        "GRENDEL_HNL_VECTORS_DIR": str(vec),
        "GRENDEL_HNL_ANALYSIS_DIR": str(analysis),
        "GRENDEL_HNL_TEMPLATES_DIR": str(templates_dir),
        "GRENDEL_HNL_CACHE_DIR": str(campaign_dir / "cache"),   # tau pool / MadGraph, shared
        "GRENDEL_FONLL_BOTTOM_GRID": str(grid_path(grid_dir, v["bottom"]["variation_tag"], "bottom")),
        "GRENDEL_FONLL_CHARM_GRID": str(grid_path(grid_dir, v["charm"]["variation_tag"], "charm")),
    })
    return env


def _run(cmd: list[str], env: dict, label: str) -> None:
    print(f"  $ {label}: {' '.join(cmd)}", flush=True)
    proc = subprocess.run(cmd, env=env, cwd=str(repo_root()))
    if proc.returncode != 0:
        raise RuntimeError(f"{label} failed (exit {proc.returncode})")


def reuse_central_channels(central_vec: Path, var_vec: Path, flavors, masses, strict: bool) -> int:
    """Hardlink the FONLL-independent channel CSVs from the central run into a
    variation. Zero-byte closed-channel sentinels are linked too. With
    ``strict`` a missing central file is an error, else skipped with a note."""
    n = 0
    for flavor in flavors:
        for label in REUSE_LABELS:
            for mass in masses:
                rel = Path(flavor) / label / f"mN_{format_mass_for_filename(mass)}.csv"
                src, dst = central_vec / rel, var_vec / rel
                if not src.exists():
                    if strict:
                        raise FileNotFoundError(f"central channel missing: {src}")
                    continue
                dst.parent.mkdir(parents=True, exist_ok=True)
                if dst.exists() or dst.is_symlink():
                    dst.unlink()
                try:
                    os.link(src, dst)
                except OSError:
                    dst.write_bytes(src.read_bytes())   # cross-device fallback
                n += 1
    return n


def run_variation(v, grid_dir, campaign_dir, source_vec, flavors, masses, args, campaign_config) -> dict:
    tag = run_tag(v, campaign_config)
    run_dir, var_vec, analysis_dir = run_dirs(campaign_dir, tag)
    env = environment(v, grid_dir, campaign_dir, tag, args.templates_dir)
    is_central = v["axis"] == "central"
    sens_csv = analysis_dir / "sensitivity.csv"
    entry = {"name": v["name"], "axis": v["axis"], "direction": v["direction"], "run_tag": tag,
             "bottom_grid_tag": v["bottom"]["variation_tag"], "charm_grid_tag": v["charm"]["variation_tag"],
             "sensitivity_csv": str(sens_csv), "configuration": campaign_config}
    if sens_csv.exists() and not args.force:
        print(f"\n=== variation {v['name']} ({v['axis']}) -> already done, skipping ===", flush=True)
        return entry
    # A "full central" run generates the MadGraph channels itself and becomes
    # the reuse source. With --reuse-from, even central only regenerates the
    # FONLL channels and hardlinks the expensive channels from the existing run.
    full_central = is_central and not args.reuse_from and not args.fonll_channels_only
    print(f"\n=== variation {v['name']} ({v['axis']}) -> {tag} ===", flush=True)
    py = [sys.executable, "-u", "-m"]
    prod = [*py, "grendel.models.hnl.production.run_all", "--flavor", *flavors, "--n-pool", str(args.n_pool),
            "--workers", str(args.prod_workers), "--seed", str(args.seed), "--skip-combine"]
    if args.mass:
        prod += ["--masses", *[str(m) for m in masses]]
    if full_central:
        prod += ["--wz-nb-core", str(args.wz_nb_core)]
    else:
        prod += ["--channels", *FONLL_CHANNELS]
    _run(prod, env, f"production[{v['name']}]")
    if not full_central and not args.fonll_channels_only:
        n = reuse_central_channels(source_vec, var_vec, flavors, masses, strict=True)
        print(f"  reused {n} channel files (Bc/Kmeson/tau/WZ) from {source_vec}", flush=True)
    combine = [*py, "grendel.models.hnl.production.combine_channels", "--flavor", *flavors]
    if args.mass:
        combine += ["--masses", *[str(m) for m in masses]]
    if args.fonll_channels_only:
        combine += ["--allow-missing"]
    _run(combine, env, f"combine[{v['name']}]")
    analyze = [*py, "grendel.models.hnl.scan", "--flavor", *flavors, "--workers", str(args.analysis_workers),
               "--decay-samples", str(args.decay_samples), "--thresholds", *map(str, args.thresholds)]
    if args.max_hit_events is not None:
        analyze += ["--max-hit-events", str(args.max_hit_events)]
    if args.event_chunk is not None:
        analyze += ["--event-chunk", str(args.event_chunk)]
    if args.analysis_seed_salt:
        analyze += ["--seed-salt", args.analysis_seed_salt]
    if args.mass:
        analyze += ["--mass", *[str(m) for m in masses]]
    _run(analyze, env, f"analysis[{v['name']}]")
    if not args.keep_intermediates:
        shutil.rmtree(var_vec, ignore_errors=True)
        print(f"  cleaned {var_vec} (kept analysis/sensitivity.csv)", flush=True)
    return entry


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grid-dir", type=Path, default=os.environ.get("GRENDEL_FONLL_GRID_DIR"),
                    required=os.environ.get("GRENDEL_FONLL_GRID_DIR") is None,
                    help="the FONLL variation campaign (variation_manifest.json + grids)")
    ap.add_argument("--campaign-dir", type=Path, default=os.environ.get("GRENDEL_HNL_CAMPAIGN_DIR"),
                    required=os.environ.get("GRENDEL_HNL_CAMPAIGN_DIR") is None,
                    help="workspace: runs/<tag>/ per variation, cache/ shared")
    ap.add_argument("--templates-dir", type=Path, default=os.environ.get("GRENDEL_HNL_TEMPLATES_DIR"),
                    required=os.environ.get("GRENDEL_HNL_TEMPLATES_DIR") is None,
                    help="decay templates shared by every variation")
    ap.add_argument("--flavor", nargs="+", default=["Ue", "Umu", "Utau"])
    ap.add_argument("--mass", type=float, nargs="+", default=None)
    ap.add_argument("--variations", nargs="+", default=None,
                    help="restrict to variation names or axes (central/scale/pdf/mass)")
    ap.add_argument("--n-pool", type=int, default=100_000)
    ap.add_argument("--seed", type=int, default=42, help="base production seed")
    ap.add_argument("--prod-workers", type=int, default=6)
    ap.add_argument("--analysis-workers", type=int, default=3)
    ap.add_argument("--decay-samples", type=int, default=100)
    ap.add_argument("--max-hit-events", type=int, default=None)
    ap.add_argument("--event-chunk", type=int, default=None)
    ap.add_argument("--analysis-seed-salt", default="")
    ap.add_argument("--thresholds", nargs="+", type=float, default=[3.0])
    ap.add_argument("--wz-nb-core", type=int, default=1)
    ap.add_argument("--reuse-from", default=None,
                    help="run tag of an existing complete run whose Bc/Kmeson/tau/WZ channels are "
                         "hardlinked into every variation, so MadGraph is never re-run")
    ap.add_argument("--fonll-channels-only", action="store_true",
                    help="plumbing test: only the four FONLL channels (incomplete physics)")
    ap.add_argument("--keep-intermediates", action="store_true",
                    help="keep each variation's llp_4vectors (tens of GB each)")
    ap.add_argument("--force", action="store_true", help="recompute finished variations")
    ap.add_argument("--registry-out", type=Path, default=None,
                    help="registry path (default: <campaign>/runs/band_registry.json)")
    args = ap.parse_args(argv)
    for name, value in (("n_pool", args.n_pool), ("decay_samples", args.decay_samples)):
        if value < 1:
            ap.error(f"--{name.replace('_', '-')} must be >= 1")

    grid_dir = args.grid_dir.expanduser().resolve()
    campaign_dir = args.campaign_dir.expanduser().resolve()
    if not (grid_dir / "variation_manifest.json").exists():
        print(f"no variation_manifest.json under {grid_dir}; produce the FONLL campaign first")
        return 1
    masses = args.mass if args.mass else list(MASS_GRID)
    variations = discover_variations(grid_dir, args.variations)
    if not variations or variations[0]["axis"] != "central":
        print("the variation set must include 'central' (the baseline to reuse from)")
        return 1
    config = campaign_configuration(args, args.flavor, masses)
    print(f"FONLL-variation band: {len(variations)} coherent runs ({', '.join(v['name'] for v in variations)})")
    print(f"  flavors={args.flavor}  masses={len(masses)}  grids={grid_dir}")
    if args.reuse_from:
        _, source_vec, _ = run_dirs(campaign_dir, args.reuse_from)
        if not source_vec.exists():
            print(f"--reuse-from run not found: {source_vec}")
            return 1
    else:
        _, source_vec, _ = run_dirs(campaign_dir, run_tag(variations[0], config))
    registry = []
    reg_path = args.registry_out or campaign_dir / "runs" / "band_registry.json"
    reg_path.parent.mkdir(parents=True, exist_ok=True)

    def write_registry():
        payload = json.dumps({"grid_dir": str(grid_dir), "configuration": config, "flavors": args.flavor,
                              "n_masses": len(masses), "variations": registry}, indent=2) + "\n"
        tmp = reg_path.with_name(f".{reg_path.name}.tmp")
        tmp.write_text(payload)
        tmp.replace(reg_path)

    for v in variations:
        registry.append(run_variation(v, grid_dir, campaign_dir, source_vec, args.flavor, masses, args, config))
        write_registry()
    print(f"\nwrote {reg_path} ({len(registry)} variations)")
    print("next: python -m grendel.band.quadrature --registry", reg_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
