"""Compare regenerated four-vector CSVs with the published ones.

    python validate_vectors.py bytes OURS_DIR REF_DIR REPORT.tsv       # files readable here (BC5)
    python validate_vectors.py sizes OURS_DIR MANIFEST PREFIX REPORT.tsv # only sizes known (BC4: inputs_manifest.txt)

``bytes``: SHA-256 equality; where the files differ, the differing lines (GNU diff) are parsed and the largest
relative difference per column reported. The CSVs print 9 significant digits (%.8e), so libm last-bit
differences between macOS arm64 and Linux x86_64 show up as +-1 in the last digit of a few lines; tolerance:
same line count, every differing line differs by < 1e-7 (weight relative to itself, E/px/py/pz relative to E,
the integer columns exactly), and < 1e-4 of the lines differ.
``sizes``: byte size equality with the manifest (row count and the sign pattern of every number).
"""
import hashlib
import subprocess
import sys
from pathlib import Path

import numpy as np


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def nlines(p):
    return int(subprocess.run(["wc", "-l", str(p)], capture_output=True, text=True).stdout.split()[0])


def by_bytes(ours, ref, report):
    rows, bad = [], []
    for f in sorted(ref.glob("mS_*.csv")):
        g = ours / f.name
        if not g.exists():
            bad.append(f"{f.name}: missing")
            continue
        same = f.stat().st_size == g.stat().st_size and sha(f) == sha(g)
        n_diff, rel, n = 0, 0.0, None
        if not same:
            n, nr = nlines(g), nlines(f)
            out = subprocess.run(["diff", str(g), str(f)], capture_output=True, text=True).stdout.splitlines()
            a = [l[2:] for l in out if l.startswith("< ")]
            b = [l[2:] for l in out if l.startswith("> ")]
            n_diff = len(a)
            if n != nr or len(a) != len(b):
                bad.append(f"{f.name}: line structure differs ({n} vs {nr} lines, {len(a)}/{len(b)} changed)")
                rel = np.inf
            elif a:
                A = np.array([[float(x) for x in l.split(",")] for l in a])
                B = np.array([[float(x) for x in l.split(",")] for l in b])
                # weight (and any extra column) relative to itself, the four-momentum relative to the energy
                d = np.abs(A - B)
                r_w = d[:, 0] / np.maximum(np.abs(B[:, 0]), 1e-300)
                r_p = d[:, 1:5].max(axis=1) / np.abs(B[:, 1])
                r_x = d[:, 5:].max(axis=1) if B.shape[1] > 5 else np.zeros(len(B))
                rel = float(np.max(np.maximum(np.maximum(r_w, r_p), r_x)))
            if rel >= 1e-7 or n_diff > 1e-4 * n:
                bad.append(f"{f.name}: {n_diff} lines differ, max rel {rel:.2g}")
        rows.append((f.name, same, n_diff, rel))
    report.write_text("file\tsha256_equal\tlines_differing\tmax_rel_diff\n" +
                      "".join(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.3g}\n" for r in rows))
    nd = np.array([r[2] for r in rows])
    print(f"{len(rows)} files: {sum(r[1] for r in rows)} byte-identical; differing lines per file: median "
          f"{np.median(nd):.0f}, max {nd.max()} (of ~1.7e6); max relative difference {max(r[3] for r in rows):.2g}")
    return bad


def by_sizes(ours, manifest, prefix, report):
    ref = {}
    for line in manifest.read_text().splitlines():
        p, s = line.split()
        if p.startswith(prefix):
            ref[Path(p).name] = int(s)
    rows, bad = [], []
    for name, size in sorted(ref.items()):
        g = ours / name
        s = g.stat().st_size if g.exists() else -1
        rows.append((name, size, s, s == size))
        if s != size:
            bad.append(f"{name}: size {s} != {size}")
    report.write_text("file\tsize_published\tsize_ours\tequal\n" + "".join(f"{a}\t{b}\t{c}\t{d}\n" for a, b, c, d in rows))
    print(f"{len(rows)} files: size equal in {sum(r[3] for r in rows)}")
    return bad


mode = sys.argv[1]
if mode == "bytes":
    bad = by_bytes(Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]))
else:
    bad = by_sizes(Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4], Path(sys.argv[5]))
for line in bad:
    print("FAIL", line)
sys.exit(1 if bad else 0)
