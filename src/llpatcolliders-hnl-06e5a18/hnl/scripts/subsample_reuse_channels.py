#!/usr/bin/env python3
"""Create a disk-light copy of a run's FONLL-independent channels for band reuse.

The full-statistics combined CSVs (~47 MB/mass at n_pool=100k) make a 21-variation
band run overflow a constrained disk. The band only measures *relative* boundary
shifts, so reduced statistics is adequate. This subsamples each reused channel
(Bc/Kmeson/tau/WZ) to at most K rows and rescales the weight column by
``n_rows / n_kept`` so the per-channel weight sum (hence the normalization and the
FONLL/independent dilution) is preserved -- only the Monte-Carlo variance grows.

Zero-byte closed-channel sentinels are copied verbatim.

    python scripts/subsample_reuse_channels.py --src full_20260606_all \
        --dst full_20260606_light --rows 10000
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np

HNL_ROOT = Path(__file__).resolve().parent.parent
if str(HNL_ROOT) not in sys.path:
    sys.path.insert(0, str(HNL_ROOT))

REUSE_LABELS = ["Bc", "Kmeson", "tau", "WZ"]
FLAVORS = ["Ue", "Umu", "Utau"]


def _runs_dir() -> Path:
    import os
    return Path(os.environ.get("HNL_TMP_DIR", HNL_ROOT / "tmp")).expanduser() / "runs"


def subsample_file(src: Path, dst: Path, k: int, rng) -> tuple[int, int]:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.stat().st_size == 0:
        dst.write_bytes(b"")          # closed-channel sentinel
        return 0, 0
    data = np.loadtxt(src, delimiter=",")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    n = len(data)
    if n <= k:
        shutil.copy2(src, dst)
        return n, n
    idx = rng.choice(n, size=k, replace=False)
    kept = data[idx].copy()
    kept[:, 0] *= n / k               # rescale weight to preserve the sum
    np.savetxt(dst, kept, delimiter=",", fmt="%.8e")
    return n, k


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", required=True, help="source run tag (e.g. full_20260606_all)")
    ap.add_argument("--dst", required=True, help="destination run tag (the light copy)")
    ap.add_argument("--rows", type=int, default=10000, help="max rows kept per channel CSV")
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args(argv)

    runs = _runs_dir()
    src_vec = runs / args.src / "llp_4vectors"
    dst_vec = runs / args.dst / "llp_4vectors"
    if not src_vec.exists():
        print(f"source not found: {src_vec}")
        return 1

    rng = np.random.default_rng(args.seed)
    n_files = tot_in = tot_out = 0
    for flavor in FLAVORS:
        for label in REUSE_LABELS:
            src_dir = src_vec / flavor / label
            if not src_dir.exists():
                continue
            for src in sorted(src_dir.glob("mN_*.csv")):
                n_in, n_out = subsample_file(src, dst_vec / flavor / label / src.name,
                                             args.rows, rng)
                n_files += 1
                tot_in += n_in
                tot_out += n_out
    print(f"subsampled {n_files} files: {tot_in:,} -> {tot_out:,} rows "
          f"(<= {args.rows}/file) into {dst_vec}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
