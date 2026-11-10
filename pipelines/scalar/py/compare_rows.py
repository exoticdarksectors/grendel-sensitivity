"""Compare scan rows (one folder per mass, or a collected table) with the published v2.5 table.

    python compare_rows.py <bc4|bc5> PUBLISHED.csv OUT_DIR_OR_TABLE [--tsv REPORT.tsv]

Per mass: n_events, n_hits, ctau, both edges and peak N at N >= 3 and N >= 10, sensitivity flags; ratios ours /
published, and whether the row is identical to the last bit. Summary per mass range: median and range of the edge
ratios. Read with the README: identical inputs give ratios 1 to ~1e-8 (the four-vectors differ from the published
ones by libm last-digit noise); the Linux exHad templates give ratios scattered by the decay Monte Carlo.
"""
import argparse
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


def show(p):
    """A path as the reports print it: relative to the repository, $GRENDEL_REF/... for the reference data."""
    p = Path(p).resolve()
    for root, name in ((os.environ.get("GRENDEL_REF", ""), "$GRENDEL_REF"),
                       (os.environ.get("REPRO") or Path(__file__).resolve().parents[3], None)):
        if root and str(root) != "none" and p.is_relative_to(Path(root).resolve()):
            rel = p.relative_to(Path(root).resolve())
            return f"{name}/{rel}" if name else str(rel)
    return str(p)

warnings.filterwarnings("ignore", message=".*float_precision.*")

ap = argparse.ArgumentParser()
ap.add_argument("bench")
ap.add_argument("published")
ap.add_argument("ours")
ap.add_argument("--tsv")
a = ap.parse_args()
rd = lambda p: pd.read_csv(p, float_precision="round_trip")  # (pandas 3 warns; still exact)
pub = rd(a.published).set_index("mass_GeV")
src = Path(a.ours)
if src.is_dir():
    parts = [rd(p) for p in sorted(src.glob("*/sensitivity.csv"))]
    ours = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=pub.reset_index().columns)
else:
    ours = rd(src)
ours = ours.set_index("mass_GeV").sort_index()
COLS = ["n_events", "n_hits", "ctau_sin2th1_m", "u2_min", "u2_max", "peak_N", "u2_min_N10", "u2_max_N10"]
rows = []
for m, r in ours.iterrows():
    if m not in pub.index:
        continue
    p = pub.loc[m]
    rec = {"mass_GeV": m}
    for c in COLS:
        x, y = float(r[c]), float(p[c])
        rec[c] = x / y if (np.isfinite(x) and np.isfinite(y) and y != 0) else (1.0 if (np.isnan(x) and np.isnan(y)) else np.nan)
    rec["sens_N3"] = f"{bool(r.has_sensitivity)}/{bool(p.has_sensitivity)}"
    rec["sens_N10"] = f"{bool(r.has_sensitivity_N10)}/{bool(p.has_sensitivity_N10)}"
    common = [c for c in pub.columns if c in ours.columns]
    rec["bit_identical"] = all((pd.isna(r[c]) and pd.isna(p[c])) or r[c] == p[c] for c in common)
    EDGES = [c for c in common if c.startswith(("u2_min", "u2_max")) and not c.endswith(
        ("_ess", "_fraction", "_open", "_open_N10"))]
    rec["edges_bit_identical"] = all((pd.isna(r[c]) and pd.isna(p[c])) or r[c] == p[c] for c in EDGES)
    rel = [abs(float(r[c]) / float(p[c]) - 1) for c in common
           if c not in ("decay_backend",) and np.isfinite(float(r[c])) and np.isfinite(float(p[c])) and float(p[c]) != 0]
    rec["max_rel_diff_all_columns"] = max(rel) if rel else 0.0
    rec["peak_N_ours"], rec["peak_N_pub"] = float(r.peak_N), float(p.peak_N)
    rows.append(rec)
t = pd.DataFrame(rows)
pd.set_option("display.width", 250)
pd.set_option("display.max_rows", 300)
fmt = {c: "{:.6f}".format for c in COLS}
print(f"{a.bench}: {len(t)} masses compared with {show(a.published)}")
if len(t):
    print(t.to_string(index=False, formatters=fmt))
    flips = t[(t.sens_N3.str.split("/").str[0] != t.sens_N3.str.split("/").str[1]) |
              (t.sens_N10.str.split("/").str[0] != t.sens_N10.str.split("/").str[1])]
    print(f"bit-identical rows: {int(t.bit_identical.sum())}/{len(t)}; rows with every edge column (u2_min*, "
          f"u2_max*) bit-identical: {int(t.edges_bit_identical.sum())}/{len(t)}; largest relative difference of any "
          f"numeric column: median {t.max_rel_diff_all_columns.median():.2g}, max {t.max_rel_diff_all_columns.max():.2g}")
    print(f"sensitivity flips (N>=3 or N>=10): "
          f"{len(flips)}" + (f" at {list(flips.mass_GeV)}" if len(flips) else ""))
    for c in ("u2_min", "u2_max", "peak_N"):
        v = t[c].dropna()
        if len(v):
            print(f"  {c}: ratio median {v.median():.4f}, range {v.min():.4f}-{v.max():.4f}, "
                  f"max |log10| {np.max(np.abs(np.log10(v))):.4f} dex")
if a.tsv:
    t.to_csv(a.tsv, sep="\t", index=False)
