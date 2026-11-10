"""Column-by-column comparison of a rebuilt table with the published one (complements compare_rows.py's per-row view).

    python compare_tables.py PUBLISHED.csv OURS.csv [--tsv REPORT.tsv] [--split MASS ...]

Works on sensitivity.csv and sensitivity_nsig10.csv alike. Rows are matched by mass_GeV. Per column: how many rows
are identical to the last bit (NaN = NaN), the largest relative difference |ours/published - 1| over the rows where
both are finite and published != 0, the mass where it occurs, and the rows where only one side is NaN. Then the
verdict columns (every *has_sensitivity* column): the masses where ours and published disagree. The file SHA-256
is compared first (bit-identical tables need no further reading). --split: masses whose largest relative difference
is reported apart from the others (09_collect.sh passes those whose decay templates differ from the published ones).
"""
import argparse
import hashlib
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
ap.add_argument("published", type=Path)
ap.add_argument("ours", type=Path)
ap.add_argument("--tsv", type=Path)
ap.add_argument("--split", nargs="*", type=float, default=[])
a = ap.parse_args()
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
rd = lambda p: pd.read_csv(p, float_precision="round_trip")  # exact parse (pandas 3 warns; still exact)
pub, ours = rd(a.published).set_index("mass_GeV"), rd(a.ours).set_index("mass_GeV")
print(f"{show(a.ours)} vs {show(a.published)}")
print(f"  SHA-256 equal: {sha(a.published) == sha(a.ours)}; rows {len(ours)} (ours) / {len(pub)} (published); "
      f"columns {len(ours.columns)} / {len(pub.columns)}")
missing = sorted(set(pub.index) - set(ours.index))
extra = sorted(set(ours.index) - set(pub.index))
if missing:
    print(f"  masses missing here: {' '.join(f'{m:g}' for m in missing)}")
if extra:
    print(f"  masses not in the published table: {' '.join(f'{m:g}' for m in extra)}")
if list(ours.columns) != list(pub.columns):
    print(f"  column sets differ: only published {sorted(set(pub.columns) - set(ours.columns))}, "
          f"only ours {sorted(set(ours.columns) - set(pub.columns))}")
idx = sorted(set(pub.index) & set(ours.index))
P, O = pub.loc[idx], ours.loc[idx]
rows = []
for c in [c for c in pub.columns if c in ours.columns]:
    p, o = P[c], O[c]
    same = (p == o) | (p.isna() & o.isna())
    rec = {"column": c, "n_rows": len(idx), "n_bit_identical": int(same.sum()), "max_rel_diff": 0.0,
           "at_mass": np.nan, "n_nan_mismatch": int((p.isna() ^ o.isna()).sum())}
    if pd.api.types.is_numeric_dtype(p) and pd.api.types.is_numeric_dtype(o) and not pd.api.types.is_bool_dtype(p):
        pf, of = p.astype(float), o.astype(float)
        ok = np.isfinite(pf) & np.isfinite(of) & (pf != 0)
        if ok.any():
            rel = (of[ok] / pf[ok] - 1).abs()
            rec["max_rel_diff"], rec["at_mass"] = float(rel.max()), float(rel.idxmax())
    rows.append(rec)
t = pd.DataFrame(rows)
ident = t[t.n_bit_identical == t.n_rows]
diff = t[t.n_bit_identical < t.n_rows].sort_values("max_rel_diff", ascending=False)
print(f"  columns bit-identical in all {len(idx)} rows: {len(ident)}/{len(t)}")
print(f"    {' '.join(ident.column)}")
print(f"  columns with differences: {len(diff)}/{len(t)} (rows differing; largest |ours/published - 1| at mass)")
for _, r in diff.iterrows():
    nan = f"; NaN on one side only in {r.n_nan_mismatch} rows" if r.n_nan_mismatch else ""
    print(f"    {r.column:34s} {r.n_rows - r.n_bit_identical:3d} rows  max {r.max_rel_diff:.3g} at {r.at_mass:g} GeV{nan}")
rowsame = ((P == O) | (P.isna() & O.isna())).all(axis=1)
print(f"  rows bit-identical in every column: {int(rowsame.sum())}/{len(idx)}")
num = [c for c in P.columns if c in O.columns and pd.api.types.is_numeric_dtype(P[c])
       and pd.api.types.is_numeric_dtype(O[c]) and not pd.api.types.is_bool_dtype(P[c])]
Pf, Of = P[num].astype(float), O[num].astype(float)
R = (Of / Pf - 1).abs().where(np.isfinite(Pf) & np.isfinite(Of) & (Pf != 0)).max(axis=1)
split = [m for m in idx if any(abs(m - x) < 1e-9 for x in a.split)]
for label, rows_ in ((f"the {len(idx) - len(split)} other masses", [m for m in idx if m not in split]),
                     (f"the split masses ({' '.join(f'{m:g}' for m in split)} GeV)", split)):
    if a.split and rows_:
        r = R.loc[rows_]
        print(f"  largest |ours/published - 1| of any numeric column, {label}: {r.max():.3g} (at {r.idxmax():g} GeV)")
verdict = [c for c in pub.columns if "has_sensitivity" in c and c in ours.columns]
flips = []
for c in verdict:
    p, o = P[c].astype(str).str.lower(), O[c].astype(str).str.lower()
    for m in p.index[p != o]:
        flips.append(f"{c} at {m:g} GeV: ours {O.at[m, c]}, published {P.at[m, c]}")
print(f"  verdict disagreements: {len(flips)} (in {len(verdict)} columns: {' '.join(verdict)})")
for f in flips:
    print(f"    {f}")
for c in ("has_sensitivity", "has_sensitivity_N10"):
    if c in O.columns:
        print(f"  {c}: {int(O[c].astype(str).str.lower().eq('true').sum())} sensitive masses here, "
              f"{int(P[c].astype(str).str.lower().eq('true').sum())} published")
if a.tsv:
    t.to_csv(a.tsv, sep="\t", index=False)
