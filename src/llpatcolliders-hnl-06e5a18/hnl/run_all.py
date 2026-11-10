#!/usr/bin/env python3
"""Parallel driver for HNL production and channel combination."""

from __future__ import annotations

import argparse
import os
import random
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

HNL_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(HNL_ROOT))
sys.path.insert(0, str(HNL_ROOT / "vendored" / "HNLCalc"))

from config_mass_grid import MASS_GRID

def _worker_meson(flavor: str, channel: str, masses: list, n_pool: int, seed: int) -> str:
    from production.decay_engine.generate_meson_csvs import (
        generate_pool, generate_bc_pool, process_channel,
    )
    rng = np.random.default_rng(seed)
    random.seed(seed)
    if channel == "bc":
        pool = generate_bc_pool(n_pool, rng)
        sigma = 0.0
    else:
        pool, sigma = generate_pool(channel, n_pool, rng)
    process_channel(flavor, channel, pool, sigma, masses, rng)
    return f"meson {flavor:>5s} {channel}"


def _worker_induced_tau(flavor: str, masses: list, n_pool: int, seed: int) -> str:
    from production.decay_engine.generate_induced_tau import build_tau_pool, process_flavor
    rng = np.random.default_rng(seed)
    random.seed(seed)
    tau_4v, tau_w, tau_asym = build_tau_pool(n_pool, rng)
    process_flavor(flavor, tau_4v, tau_w, tau_asym, masses, rng)
    return f"itau  {flavor:>5s}"


def _worker_kaon(flavor: str, masses: list, n_pool: int, seed: int) -> str:
    from production.decay_engine.generate_kaon_csvs import sample_kaon_4vectors, process_flavor
    rng = np.random.default_rng(seed)
    random.seed(seed)
    pool = sample_kaon_4vectors(n_pool, rng)
    process_flavor(flavor, pool, masses, rng)
    return f"kaon  {flavor:>5s}"


def _worker_baryon(flavor, masses, n_pool, seed):
    from production.decay_engine.generate_baryon_csvs import build_lambda_b_pool, process_flavor
    rng = np.random.default_rng(seed); random.seed(seed)
    pool_v, sigma_b = build_lambda_b_pool(n_pool, rng)
    process_flavor(flavor, pool_v, sigma_b, masses, rng)
    return f"baryon {flavor:>5s}"


def _worker_prompt_tau_stage2(flavor: str, masses: list, seed: int) -> str:
    from production.madgraph.run_tau_production import _load_pool, process_flavor
    rng = np.random.default_rng(seed)
    random.seed(seed)
    pool_w, pool_E, pool_px, pool_py, pool_pz, pool_origin = _load_pool()
    process_flavor(flavor, pool_w, pool_E, pool_px, pool_py, pool_pz, pool_origin,
                   masses, rng)
    return f"ptau  {flavor:>5s}"


def _worker_wz(flavor: str, masses: list, n_events: int, nb_core: int) -> str:
    from production.madgraph.run_wz_production import run_single_point
    failed = []
    for m in masses:
        if not run_single_point(flavor, m, n_events, nb_core=nb_core):
            failed.append(m)
    if failed:
        raise RuntimeError(f"W/Z failed for {flavor} masses: {failed}")
    return f"wz    {flavor:>5s} ({len(masses)}/{len(masses)} points)"


def main():
    ap = argparse.ArgumentParser(description="Parallel HNL production driver")
    ap.add_argument("--flavor", choices=["Ue", "Umu", "Utau"], nargs="+",
                    default=["Ue", "Umu", "Utau"])
    ap.add_argument("--n-pool", type=int, default=100_000,
                    help="Sampled mesons/taus per (flavor, channel) job (default: 100k)")
    ap.add_argument("--workers", type=int, default=None,
                    help="Max concurrent processes (default: os.cpu_count())")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--masses", type=float, nargs="+", default=None,
                    help="Custom mass list (default: full grid)")
    ap.add_argument("--skip-combine", action="store_true",
                    help="Skip the final combine_channels step")
    ap.add_argument("--no-wz", action="store_true",
                    help="Skip the W/Z -> l N MadGraph channel (default on)")
    ap.add_argument("--wz-nevents", type=int, default=None,
                    help="Events per (flavor, mass) MadGraph point "
                         "(default: config N_EVENTS_DEFAULT)")
    ap.add_argument("--wz-nb-core", type=int, default=1,
                    help="CPU cores per MadGraph generate_events job")
    ap.add_argument("--no-prompt-tau", action="store_true",
                    help="Skip the prompt-tau channel (Stage 1 MG5 + Stage 2 decay)")
    ap.add_argument("--prompt-tau-nevents", type=int, default=None,
                    help="Stage 1 tau-pool size (default: config N_EVENTS_DEFAULT). "
                         "Ignored if a compatible cached tau pool already exists.")
    ap.add_argument("--prompt-tau-nb-core", type=int, default=1,
                    help="CPU cores for Stage 1 MG5 generate_events")
    ap.add_argument("--channels", nargs="+", default=None,
                    choices=["bottom", "charm", "bc", "induced_tau",
                             "kaon", "baryon", "prompt_tau", "wz"],
                    help="Generate only these channels (default: all). Used by the "
                         "FONLL-variation band driver to regenerate only the "
                         "FONLL-dependent channels.")
    args = ap.parse_args()

    workers = args.workers or os.cpu_count() or 1
    masses = args.masses if args.masses else MASS_GRID

    # Channel selection: --channels narrows the set; --no-wz / --no-prompt-tau
    # still subtract from whatever was selected (backward compatible).
    all_channels = ["bottom", "charm", "bc", "induced_tau",
                    "kaon", "baryon", "prompt_tau", "wz"]
    selected = set(args.channels) if args.channels else set(all_channels)
    if args.no_wz:
        selected.discard("wz")
    if args.no_prompt_tau:
        selected.discard("prompt_tau")

    meson_channels = [ch for ch in ("bottom", "charm", "bc") if ch in selected]
    from production.paths import describe_paths
    print("HNL artifact paths:")
    print(describe_paths())

    if "prompt_tau" in selected:
        from production.madgraph.run_tau_production import POOL_CSV, generate_tau_pool
        from production.paths import existing_tau_pool_csv
        from config_mass_grid import N_EVENTS_DEFAULT
        pt_nevents = args.prompt_tau_nevents or N_EVENTS_DEFAULT

        regen_reason = None
        pool_csv = existing_tau_pool_csv()
        if not pool_csv.exists() or pool_csv.stat().st_size == 0:
            regen_reason = "no cached pool"
        else:
            with open(pool_csv) as f:
                n_rows = sum(1 for _ in f)
            threshold = max(int(0.5 * pt_nevents), 1)
            if n_rows < threshold:
                regen_reason = (f"cached pool has {n_rows} rows < threshold "
                                f"{threshold} for --prompt-tau-nevents={pt_nevents}")
            else:
                print(f"Prompt-tau Stage 1: reusing cached pool at {pool_csv} "
                      f"({n_rows} rows)")

        if regen_reason is not None:
            print(f"Prompt-tau Stage 1: regenerating ({regen_reason}); "
                  f"building pool ({pt_nevents} events) at {POOL_CSV} ...")
            ok = generate_tau_pool(pt_nevents, nb_core=args.prompt_tau_nb_core)
            if not ok:
                print("Prompt-tau Stage 1 failed; aborting.")
                return 1

    jobs = []
    seed = args.seed
    for flavor in args.flavor:
        for ch in meson_channels:
            jobs.append((_worker_meson, (flavor, ch, masses, args.n_pool, seed)))
            seed += 1
        if "induced_tau" in selected:
            jobs.append((_worker_induced_tau, (flavor, masses, args.n_pool, seed)))
            seed += 1
        if "kaon" in selected:
            jobs.append((_worker_kaon, (flavor, masses, args.n_pool, seed)))
            seed += 1
        if "baryon" in selected:
            jobs.append((_worker_baryon, (flavor, masses, args.n_pool, seed)))
            seed += 1
        if "prompt_tau" in selected:
            jobs.append((_worker_prompt_tau_stage2, (flavor, masses, seed)))
            seed += 1

    if "wz" in selected:
        from config_mass_grid import N_EVENTS_DEFAULT
        wz_nevents = args.wz_nevents or N_EVENTS_DEFAULT
        for flavor in args.flavor:
            jobs.append((_worker_wz, (flavor, masses, wz_nevents, args.wz_nb_core)))

    print(f"Submitting {len(jobs)} jobs across {workers} workers "
          f"(flavors={args.flavor}, masses={len(masses)}, n_pool={args.n_pool})")
    t0 = time.time()
    failures = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(fn, *a): a for fn, a in jobs}
        for f in as_completed(futures):
            try:
                label = f.result()
                print(f"  [{time.time()-t0:6.1f}s]  done: {label}")
            except Exception as exc:
                a = futures[f]
                failures.append((a, exc))
                print(f"  FAILED job args={a}: {exc}")

    if failures:
        print(f"\n{len(failures)} job(s) failed.")
        return 1

    if not args.skip_combine:
        print("\nCombining channels...")
        from production.combine_channels import combine_for_point
        active_channels = ["Bmeson", "Dmeson", "Bc", "Bbaryon", "induced_tau", "Kmeson"]
        if not args.no_prompt_tau:
            active_channels.append("tau")
        if not args.no_wz:
            active_channels.append("WZ")
        n_combined = 0
        for flavor in args.flavor:
            for m in masses:
                if combine_for_point(flavor, m, channels=active_channels,
                                     strict=True) > 0:
                    n_combined += 1
        print(f"  combined files written for {n_combined} (flavor, mass) points")

    print(f"\nTotal wall time: {time.time()-t0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
