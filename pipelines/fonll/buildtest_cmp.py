"""Compare raw quark-grid rows of buildtest_rows.py with the stored macOS raw grid (same y nodes); see buildtest_rows.py."""
import sys
from pathlib import Path
import numpy as np

test, ref = Path(sys.argv[1]), Path(sys.argv[2])
iys = [int(x) for x in sys.argv[3].split(",")]
npt = 80
t = test.read_text().splitlines()
r = ref.read_text().splitlines()
print("header equal:", t[0] == r[0])
rows_t = t[1:]
rows_r = [r[1 + i * npt + j] for i in iys for j in range(npt)]
assert len(rows_t) == len(rows_r), (len(rows_t), len(rows_r))
same = sum(a == b for a, b in zip(rows_t, rows_r))
print(f"identical lines: {same}/{len(rows_r)}")
a = np.array([[float(x) for x in l.split()] for l in rows_t])
b = np.array([[float(x) for x in l.split()] for l in rows_r])
assert np.array_equal(a[:, :2], b[:, :2])
tot_a = a[:, 2:6].sum(axis=1)   # sum of the four FONLL components (cols 3-6)
tot_b = b[:, 2:6].sum(axis=1)
nz = tot_b != 0
rel = np.abs(tot_a[nz] - tot_b[nz]) / np.abs(tot_b[nz])
print(f"sum of components: max rel diff {rel.max():.3e}, median {np.median(rel):.3e}")
for c in range(2, 6):
    m = b[:, c] != 0
    rc = np.abs(a[m, c] - b[m, c]) / np.abs(b[m, c])
    print(f"  col {c + 1}: max rel diff {rc.max():.3e}")
