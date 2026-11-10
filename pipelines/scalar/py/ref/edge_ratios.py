"""Edges of the comparison curves against the GRENDEL table, for the numbers quoted in Sec. IV.

    python edge_ratios.py bc4 <sensitivity.csv>
    python edge_ratios.py hnl <sensitivity.csv>

A projection is a closed contour: at each GRENDEL mass inside it, its lower (upper) edge is the smallest (largest)
coupling where the contour crosses that mass. Ratios are projection / GRENDEL; run with PYTHONPATH = the code.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from curves.plot.common import load_polygon_segments

D = Path(__import__("curves").__file__).resolve().parent / "data"


def crossings(segments, x):
    ys = []
    for s in segments:
        s = s[np.isfinite(s).all(axis=1)]
        for (x0, y0), (x1, y1) in zip(s[:-1], s[1:]):
            if (x0 - x) * (x1 - x) <= 0 and x0 != x1:
                t = (x - x0) / (x1 - x0)
                ys.append(10 ** (np.log10(y0) + t * (np.log10(y1) - np.log10(y0))))
    return (min(ys), max(ys)) if len(ys) >= 2 else (np.nan, np.nan)


def table(rows, segs, name, lo=None, hi=None):
    out = []
    for _, r in rows.iterrows():
        a, b = crossings(segs, r.mass_GeV)
        out.append((r.mass_GeV, a / r.u2_min, b / r.u2_max, a, b))
    t = pd.DataFrame(out, columns=["m", "lower_ratio", "upper_ratio", "proj_lower", "proj_upper"]).dropna()
    if lo is not None:
        t = t[t.m.between(lo, hi)]
    return t


def summary(t, name, ranges):
    for lo, hi in ranges:
        s = t[t.m.between(lo, hi)]
        if len(s):
            print(f"  {name:8s} {lo:>5}-{hi:<5} GeV ({len(s):2d}): lower proj/GRENDEL {s.lower_ratio.min():.2f}-"
                  f"{s.lower_ratio.max():.2f} (median {s.lower_ratio.median():.2f}); upper {s.upper_ratio.min():.2f}-"
                  f"{s.upper_ratio.max():.2f} (median {s.upper_ratio.median():.2f})")


bench, csv = sys.argv[1], sys.argv[2]
g = pd.read_csv(csv, float_precision="round_trip")
if bench == "bc4":
    g = g[g.has_sensitivity].sort_values("mass_GeV")
    print(f"BC4: deepest {g.u2_min.min():.3g} at {g.loc[g.u2_min.idxmin(), 'mass_GeV']} GeV; reach N>=3 "
          f"{g.mass_GeV.max()} GeV; lower edge at 1-3.5 GeV {g[g.mass_GeV.between(1, 3.5)].u2_min.min():.3g}-"
          f"{g[g.mass_GeV.between(1, 3.5)].u2_min.max():.3g}")
    n10 = g[g.has_sensitivity_N10] if "has_sensitivity_N10" in g else g.iloc[:0]
    print(f"     reach N>=10 {n10.mass_GeV.max() if len(n10) else 'n/a'} GeV")
    P = D / "bc4" / "projections"
    for name, f in (("FASER2", "faser2_bc4.dat"), ("SHiP", "ship_bb_bc4.dat"), ("CODEX-b", "codexb_300fb.dat")):
        summary(table(g, load_polygon_segments(P / f), name), name, ((0.1, 1.0), (1.0, 1.5), (1.5, 3.0), (3.0, 4.0), (0.1, 5)))
elif bench == "hnl":
    for fl in ("Ue", "Umu", "Utau"):
        s = g[(g.flavor == fl) & g.has_sensitivity].sort_values("mass_GeV")
        k = s.u2_min.idxmin()
        print(f"{fl}: deepest {s.u2_min.min():.3g} at {s.loc[k, 'mass_GeV']} GeV; at 1.5 GeV "
              f"{s[np.isclose(s.mass_GeV, 1.5)].u2_min.min():.3g}; at 0.4 GeV {s[np.isclose(s.mass_GeV, 0.4)].u2_min.min():.3g}; "
              f"closure (last sensitive mass) {s.mass_GeV.max()} GeV")
else:
    sys.exit("bc4 or hnl")
