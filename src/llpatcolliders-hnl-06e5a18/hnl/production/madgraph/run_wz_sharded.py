#!/usr/bin/env python3
"""Run W/Z HNL production in mass shards with isolated MG5 work dirs."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config_mass_grid import MASS_GRID, N_EVENTS_DEFAULT
from production.io import llp_csv_path
from production.paths import CACHE_DIR, TMP_DIR

RUN_WZ = Path(__file__).with_name("run_wz_production.py")


def _mass_label(mass: float) -> str:
    return f"{mass:.3f}".rstrip("0").rstrip(".")


def _round_robin_chunks(items: list[float], n_chunks: int) -> list[list[float]]:
    chunks = [[] for _ in range(n_chunks)]
    for i, item in enumerate(items):
        chunks[i % n_chunks].append(item)
    return [chunk for chunk in chunks if chunk]


def _run_shard(payload: tuple) -> tuple[int, int, Path, float]:
    shard_id, flavor, masses, nevents, nb_core, work_base, log_dir = payload
    shard_work = work_base / f"shard_{shard_id:02d}"
    log_path = log_dir / f"wz_{flavor}_shard_{shard_id:02d}.log"
    shard_work.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable, "-u", str(RUN_WZ),
        "--flavor", flavor,
        "--nevents", str(nevents),
        "--nb-core", str(nb_core),
        "--masses",
        *[_mass_label(m) for m in masses],
    ]
    env = os.environ.copy()
    env["HNL_MG5_WORK_DIR"] = str(shard_work)
    env["PYTHONUNBUFFERED"] = "1"

    t0 = time.time()
    with log_path.open("w") as log:
        log.write(f"CMD: {' '.join(cmd)}\n")
        log.write(f"HNL_MG5_WORK_DIR={shard_work}\n\n")
        log.flush()
        result = subprocess.run(
            cmd,
            cwd=PROJECT_ROOT,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    return shard_id, result.returncode, log_path, time.time() - t0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Parallel W/Z -> l N production over mass shards")
    parser.add_argument("--flavor", choices=["Ue", "Umu", "Utau"], default="Umu")
    parser.add_argument("--masses", type=float, nargs="+", default=None)
    parser.add_argument("--min-mass", type=float, default=None)
    parser.add_argument("--nevents", type=int, default=N_EVENTS_DEFAULT)
    parser.add_argument("--jobs", type=int, default=8,
                        help="Number of independent mass shards")
    parser.add_argument("--nb-core", type=int, default=1,
                        help="MG5 cores per shard; jobs * nb-core is the target core count")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip non-empty WZ CSVs that already exist")
    parser.add_argument("--work-base", type=Path,
                        default=CACHE_DIR / "madgraph_wz_shards",
                        help="Base directory for per-shard MG5 work dirs")
    parser.add_argument("--log-dir", type=Path,
                        default=TMP_DIR / "wz_shards",
                        help="Directory for per-shard logs")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    if args.jobs < 1:
        print("ERROR: --jobs must be >= 1")
        return 1
    if args.nb_core < 1:
        print("ERROR: --nb-core must be >= 1")
        return 1

    masses = args.masses if args.masses else MASS_GRID
    if args.min_mass is not None:
        masses = [m for m in masses if m >= args.min_mass]
    if args.skip_existing:
        masses = [
            m for m in masses
            if not (
                (path := llp_csv_path(args.flavor, "WZ", m, mkdir=False)).exists()
                and path.stat().st_size > 0
            )
        ]

    if not masses:
        print("No W/Z mass points to run.")
        return 0

    chunks = _round_robin_chunks(masses, min(args.jobs, len(masses)))

    print("Sharded W/Z production")
    print(f"  Flavor: {args.flavor}")
    print(f"  Masses: {len(masses)} points ({min(masses):.3f} - {max(masses):.3f} GeV)")
    print(f"  Jobs: {len(chunks)}")
    print(f"  MG5 cores/job: {args.nb_core}")
    print(f"  Target core use: {len(chunks) * args.nb_core}")
    print(f"  Work base: {args.work_base}")
    print(f"  Logs: {args.log_dir}")
    for i, chunk in enumerate(chunks):
        print(
            f"    shard {i:02d}: {len(chunk):3d} masses "
            f"({chunk[0]:.3f} ... {chunk[-1]:.3f})")

    if args.dry_run:
        return 0

    t0 = time.time()
    payloads = [
        (i, args.flavor, chunk, args.nevents, args.nb_core,
         args.work_base, args.log_dir)
        for i, chunk in enumerate(chunks)
    ]

    failures = []
    with ProcessPoolExecutor(max_workers=len(chunks)) as pool:
        futures = {pool.submit(_run_shard, payload): payload[0] for payload in payloads}
        for future in as_completed(futures):
            shard_id, code, log_path, elapsed = future.result()
            if code == 0:
                print(f"  done shard {shard_id:02d} ({elapsed:.1f}s): {log_path}")
            else:
                print(f"  FAILED shard {shard_id:02d} ({elapsed:.1f}s): {log_path}")
                failures.append((shard_id, log_path))

    print(f"\nTotal wall time: {time.time() - t0:.1f}s")
    if failures:
        print("Failed shards:")
        for shard_id, log_path in failures:
            print(f"  shard {shard_id:02d}: {log_path}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
