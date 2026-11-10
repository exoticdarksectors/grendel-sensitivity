"""Pixel comparison of two PNGs (a rendered Fig. 6 against the Overleaf file).

    python compare_png.py <a.png> <b.png> [--diff out.png]

Prints the sizes, the number of differing pixels (any RGBA channel), their bounding box, the largest channel
difference and, with --diff, writes an image of b dimmed with the differing pixels in red. Paths inside the
repository are printed relative to it.
"""
import argparse
import hashlib
import os
from pathlib import Path

import numpy as np
from PIL import Image

ap = argparse.ArgumentParser()
ap.add_argument("a", type=Path)
ap.add_argument("b", type=Path)
ap.add_argument("--diff", type=Path, default=None)
a = ap.parse_args()
x = np.asarray(Image.open(a.a).convert("RGBA")).astype(int)
y = np.asarray(Image.open(a.b).convert("RGBA")).astype(int)
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()[:16]
root = Path(os.environ.get("REPRO", "/")).resolve()
show = lambda p: p.resolve().relative_to(root) if p.resolve().is_relative_to(root) else p
print(f"a: {show(a.a)} {x.shape[1]}x{x.shape[0]} sha256 {sha(a.a)}...")
print(f"b: {show(a.b)} {y.shape[1]}x{y.shape[0]} sha256 {sha(a.b)}...")
if x.shape != y.shape:
    print("sizes differ: no pixel comparison")
    raise SystemExit(0)
d = np.abs(x - y).max(axis=-1)
n = int((d > 0).sum())
print(f"differing pixels: {n} of {d.size} ({100 * n / d.size:.3f} %); largest channel difference {int(d.max())}")
if n:
    r, c = np.nonzero(d)
    print(f"bounding box of the differences: x {c.min()}-{c.max()}, y {r.min()}-{r.max()}")
    for t in (8, 32, 128):
        print(f"  pixels differing by more than {t}: {int((d > t).sum())}")
    if a.diff:
        img = (y[..., :3] * 0.3 + 255 * 0.7).astype(np.uint8)
        img[d > 0] = [255, 0, 0]
        Image.fromarray(img).save(a.diff)
        print(f"diff image: {show(a.diff)}")
