#!/usr/bin/env python3
"""Every number Sec. IV.C (BC5) of the GRENDEL paper quotes, from one sensitivity table.

    python paper_numbers.py TABLE.csv [CURVES_ROOT]

CURVES_ROOT is a grendel-sensitivity checkout (the projections and the LHC bands); default: src/grendel-sensitivity
of this repository. (The only change from the campaign copy, bc5_v2_1_20261003/results/bc5_hybrid/paper_numbers.py,
is that default, which pointed at the original machine's reference tree.)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.path import Path as MPath

table = Path(sys.argv[1])
root = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).resolve().parents[4] / "src/grendel-sensitivity"
g = pd.read_csv(table).set_index("mass_GeV").sort_index()
g = g[g.index >= 0.22]


def chains(path):
    out, cur = [], []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line.startswith("#"):
            continue
        if not line:
            if cur:
                out.append(np.array(cur))
                cur = []
            continue
        cur.append([float(v) for v in line.split()[:2]])
    if cur:
        out.append(np.array(cur))
    return out


def edges(ch, m):
    ys = []
    for c in ch:
        x, y, lm = np.log10(c[:, 0]), np.log10(c[:, 1]), np.log10(m)
        for i in range(len(c) - 1):
            if (x[i] - lm) * (x[i + 1] - lm) <= 0 and x[i] != x[i + 1]:
                ys.append(y[i] + (y[i + 1] - y[i]) * (lm - x[i]) / (x[i + 1] - x[i]))
    return (10 ** min(ys), 10 ** max(ys)) if ys else (np.nan, np.nan)


proj = {k: chains(root / "curves/data/bc5/projections" / f) for k, f in
        [("SHiP", "ship_pbc_bc5.dat"), ("CODEX-b", "codexb_pbc_bc5.dat"),
         ("ANUBIS", "anubis_pbc_bc5.dat"), ("MATHUSLA40", "mathusla40_pbc_bc5.dat")]}
ratio = {}
for k, ch in proj.items():
    rows = []
    for m, r in g.iterrows():
        lo, hi = edges(ch, m)
        rows.append((m, lo / r.u2_min, hi / r.u2_max))
    ratio[k] = pd.DataFrame(rows, columns=["m", "lo", "hi"]).set_index("m")

print(f"table: {table}")
for m in (1.0, 10.0, 60.0):
    print(f"  island at {m:g} GeV: {g.loc[m, 'u2_min']:.2e} - {g.loc[m, 'u2_max']:.2e}")
k = g.u2_min.idxmin()
print(f"  smallest mixing: {g.u2_min.min():.2e} at {k:g} GeV")
print(f"  peak N minimum: {g.peak_N.min():.0f} at {g.peak_N.idxmin():g} GeV; N10 everywhere: {bool(g.has_sensitivity_N10.all())}")
print(f"  N>=10: lower edge x{(g.u2_min_N10 / g.u2_min).max():.2f} max; upper edge /{(g.u2_max / g.u2_max_N10).min():.2f}-"
      f"{(g.u2_max / g.u2_max_N10).max():.2f} (median {(g.u2_max / g.u2_max_N10).median():.2f})")

s = ratio["SHiP"]
lo = s.loc[0.3:2.5, "lo"]
print(f"  SHiP below 2.5 GeV: probes {1 / lo.max():.1f}-{1 / lo.min():.1f}x smaller mixing")
deeper = s[(s.index > 2) & (s.lo > 1)]
print(f"  GRENDEL deeper than SHiP from {deeper.index.min():g} GeV ({deeper.lo.min():.2f}x there) to its end")
c = ratio["CODEX-b"]
cross = c[(c.index > 2) & (c.lo > 1)].index.min()
print(f"  CODEX-b lower edge: up to {1 / c.loc[:cross, 'lo'].min():.1f}x smaller below ~{cross:g} GeV; "
      f"GRENDEL up to {c.loc[cross:, 'lo'].max():.1f}x smaller above")
up = c.loc[1.0:, "hi"]
print(f"  CODEX-b upper edge above 1 GeV: GRENDEL reaches >= {1 / up.max():.2f}x larger (at {up.idxmax():g} GeV)")
mt = ratio["MATHUSLA40"].loc[10.5:62]
print(f"  MATHUSLA40 10-60 GeV: lower edge {1 / mt.lo.max():.1f}-{1 / mt.lo.min():.1f}x smaller; "
      f"GRENDEL upper edge {1 / mt.hi.max():.1f}-{1 / mt.hi.min():.1f}x larger")
an = ratio["ANUBIS"].loc[5.5:62]
print(f"  ANUBIS >5 GeV: lower edge {1 / an.lo.max():.2f}-{1 / an.lo.min():.1f}x smaller; upper edge {an.hi.min():.2f}-{an.hi.max():.1f}x larger")

bands = []
for stem in ("lhc_hss_hadronic_2402.01898", "lhc_hss_mumu_2112.13769", "lhc_hss_mumu_2402.14491",
             "lhc_hss_bb_llp_summaries"):
    s2 = chains(root / "curves/data/bc5/existing" / f"{stem}.dat")
    for a, b in zip(s2[0::2], s2[1::2]):
        p = np.log10(np.r_[a, b[::-1], a[:1]])
        bands.append(MPath(p, closed=True))
cov = []
for m in g.loc[2.0:55.0].index:
    r = g.loc[m]
    ys = np.linspace(np.log10(r.u2_min), np.log10(r.u2_max), 400)
    pts = np.c_[np.full_like(ys, np.log10(m)), ys]
    hit = np.zeros(len(ys), bool)
    for b in bands:
        hit |= b.contains_points(pts)
    cov.append((m, hit.mean() * (ys[-1] - ys[0]), hit.mean()))
cov = np.array(cov)
print(f"  LHC h->SS bands cover the upper {cov[:, 1].min():.1f}-{cov[:, 1].max():.1f} decades "
      f"({cov[:, 2].min():.0%}-{cov[:, 2].max():.0%}) of the island at 2-55 GeV")
