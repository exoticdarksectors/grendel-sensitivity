#!/usr/bin/env python3
"""Combine per-channel HNL CSVs into per-mass files."""

import sys
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config_mass_grid import MASS_GRID
from production.io import llp_csv_path, read_csv_matrix, write_csv_matrix
from production.paths import LLP_VECTORS_DIR

OUTPUT_BASE = LLP_VECTORS_DIR

CHANNELS = ["Bmeson", "Dmeson", "Bc", "Bbaryon", "tau", "induced_tau", "Kmeson", "WZ"]
FLAVORS = ["Ue", "Umu", "Utau"]


def combine_for_point(flavor, mass, channels=None, strict=True):
    combined_path = llp_csv_path(flavor, "combined", mass, base=OUTPUT_BASE)

    channels = CHANNELS if channels is None else channels
    all_data, missing, bad = [], [], []
    for ch in channels:
        csv_path = llp_csv_path(flavor, ch, mass, base=OUTPUT_BASE, mkdir=False)
        if not csv_path.exists():
            missing.append(ch)
            continue
        if csv_path.stat().st_size == 0:
            continue
        try:
            data = read_csv_matrix(csv_path)
            if data.shape[1] >= 5:
                all_data.append(data[:, :5])
            else:
                bad.append(f"{ch}: only {data.shape[1]} columns")
        except Exception as e:
            bad.append(f"{ch}: {e}")

    if strict and (missing or bad):
        mass_label = combined_path.stem.removeprefix("mN_")
        raise RuntimeError(
            f"incomplete combine for {flavor} mN={mass_label}: "
            f"missing={missing}, bad={bad}"
        )

    if all_data:
        combined = np.vstack(all_data)
        write_csv_matrix(combined_path, combined)
        return len(combined)
    else:
        write_csv_matrix(combined_path, [])
        return 0


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Combine HNL production channels")
    parser.add_argument("--flavor", choices=FLAVORS, nargs="+", default=FLAVORS)
    parser.add_argument("--masses", type=float, nargs="+", default=None)
    parser.add_argument("--allow-missing", action="store_true",
                        help="skip missing or malformed channels instead of failing "
                             "(for deliberate partial runs; zero-byte closed-channel "
                             "sentinels are always accepted)")
    args = parser.parse_args()

    masses = args.masses if args.masses else MASS_GRID

    total_files = 0
    total_events = 0
    n_empty = 0

    for flavor in args.flavor:
        print(f"\n{flavor}:")
        for mass in masses:
            n_ev = combine_for_point(flavor, mass, strict=not args.allow_missing)
            total_files += 1
            total_events += n_ev
            if n_ev == 0:
                n_empty += 1

        combined_dir = llp_csv_path(flavor, "combined", masses[0], base=OUTPUT_BASE).parent
        n_files = len(list(combined_dir.glob("mN_*.csv")))
        print(f"  {n_files} files in {combined_dir}")

    print(f"\nTotal: {total_files} combined files, {total_events} events, {n_empty} empty")


if __name__ == "__main__":
    main()
