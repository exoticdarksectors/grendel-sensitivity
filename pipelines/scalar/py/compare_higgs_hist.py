"""Compare a rebuilt Higgs (pT, y) histogram with the tracked one.

    python compare_higgs_hist.py OURS.csv.gz TRACKED.csv.gz

Per mode: cells and counts identical?, the per-mode summary line of the header, and a chi^2 of the two count
tables (two-sample, cells merged until each holds >= 10 counts in the sum) with its p-value, plus the
projections' means. Prints a table; exit status 0 always (the caller decides how to read it).
"""
import gzip
import sys

import numpy as np
from scipy import stats


def read(path):
    head, cells = [], {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                head.append(line.strip())
                continue
            if line.strip():
                mode, a, b, c, d, n = line.strip().split(",")
                cells.setdefault(mode, {})[(float(a), float(b), float(c), float(d))] = int(n)
    return head, cells


ha, A = read(sys.argv[1])
hb, B = read(sys.argv[2])
same_lines = [x for x in ha if not x.startswith("# built")] == [x for x in hb if not x.startswith("# built")]
print(f"header identical apart from the build time: {same_lines}")
for mode in B:
    a, b = A.get(mode, {}), B[mode]
    keys = sorted(set(a) | set(b))
    na = np.array([a.get(k, 0) for k in keys], float)
    nb = np.array([b.get(k, 0) for k in keys], float)
    identical = a == b
    # merge sparse cells (ordered by pT then y) until each merged cell holds >= 10 counts in total
    groups, cur_a, cur_b = [], 0.0, 0.0
    for x, y in zip(na, nb):
        cur_a += x; cur_b += y
        if cur_a + cur_b >= 10:
            groups.append((cur_a, cur_b)); cur_a = cur_b = 0.0
    if cur_a + cur_b and groups:
        ga, gb = groups[-1]; groups[-1] = (ga + cur_a, gb + cur_b)
    g = np.array(groups)
    sa, sb = g[:, 0].sum(), g[:, 1].sum()
    k1, k2 = np.sqrt(sb / sa), np.sqrt(sa / sb)
    chi2 = float(np.sum((k1 * g[:, 0] - k2 * g[:, 1]) ** 2 / (g[:, 0] + g[:, 1])))
    ndf = len(g) - 1
    p = float(stats.chi2.sf(chi2, ndf))
    pt = np.array([(k[0] + k[1]) / 2 for k in keys])
    y = np.array([abs(k[2] + k[3]) / 2 for k in keys])
    line_a = next((x for x in ha if x.startswith(f"# {mode}:")), "")
    line_b = next((x for x in hb if x.startswith(f"# {mode}:")), "")
    print(f"{mode:4s} cells {len(a)}/{len(b)} counts {int(na.sum())}/{int(nb.sum())} identical {identical}; "
          f"chi2/ndf {chi2:.1f}/{ndf} p={p:.3f}; binned <pT> {np.average(pt, weights=na):.2f}/"
          f"{np.average(pt, weights=nb):.2f} GeV, <|y|> {np.average(y, weights=na):.3f}/{np.average(y, weights=nb):.3f}")
    if line_a != line_b:
        print(f"     ours:    {line_a[2:]}\n     tracked: {line_b[2:]}")
