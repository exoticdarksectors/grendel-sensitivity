"""Pixel comparison of two PNGs.

    python compare_png.py A.png B.png [DIFF.png]

Prints the sizes, the fraction of pixels that differ (any channel), the fraction differing by more than 32/255,
the largest channel difference and the bounding box of the differing pixels; writes a diff image (differing
pixels red on a faded copy of B) if asked.
"""
import os
import sys
from pathlib import Path

import numpy as np
from PIL import Image


def show(p):
    """A path as the reports print it: relative to the repository, $GRENDEL_REF/... for the reference data."""
    p = Path(p).resolve()
    for root, name in ((os.environ.get("GRENDEL_REF", ""), "$GRENDEL_REF"),
                       (os.environ.get("REPRO") or Path(__file__).resolve().parents[3], None)):
        if root and str(root) != "none" and p.is_relative_to(Path(root).resolve()):
            rel = p.relative_to(Path(root).resolve())
            return f"{name}/{rel}" if name else str(rel)
    return str(p)

a = np.asarray(Image.open(sys.argv[1]).convert("RGBA")).astype(int)
b = np.asarray(Image.open(sys.argv[2]).convert("RGBA")).astype(int)
print(f"{show(sys.argv[1])}: {a.shape[1]}x{a.shape[0]}; {show(sys.argv[2])}: {b.shape[1]}x{b.shape[0]}")
if a.shape != b.shape:
    print("different sizes: no pixel comparison")
    sys.exit(0)
d = np.abs(a - b).max(axis=2)
n = d.size
ys, xs = np.nonzero(d)
box = f"x {xs.min()}-{xs.max()}, y {ys.min()}-{ys.max()}" if len(xs) else "none"
print(f"pixels differing: {np.count_nonzero(d) / n:.5%} (>32/255: {np.count_nonzero(d > 32) / n:.5%}); "
      f"max channel difference {d.max()}; bounding box {box}")
if len(sys.argv) > 3:
    base = (255 - (255 - b[..., :3]) * 0.3).astype(np.uint8)
    base[d > 0] = [255, 0, 0]
    Image.fromarray(base).save(sys.argv[3])
