"""Row-by-row comparison of a rebuilt HNL table with the published one (curves/data/grendel_v2_3/hnl).

    python compare_table.py <rebuilt sensitivity.csv> <published sensitivity.csv> [--out DIR] [--tag NAME]

Rows are matched on (flavor, mass_GeV); block = coarse (grids.COARSE, the lxplus rows) or fine (grids.FINE, the
HTCondor rows). Reports
  * per column: rows whose CSV text is identical (bit-identical value) and rows whose parsed value is identical
    (NaN == NaN), overall and per block;
  * per numeric column, flavour and block: the maximum and median |new/published - 1| over rows where both are
    finite, and the rows where only one side is empty (an edge that closes in one table and not in the other);
  * every row whose verdict changes: has_sensitivity, has_sensitivity_N10 and the open-edge flags.
Writes <DIR>/table_<tag>_rows.csv (per-row relative differences) and <DIR>/table_<tag>.md, and prints the markdown.
"""
import argparse
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import grids  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("new", type=Path)
ap.add_argument("ref", type=Path)
ap.add_argument("--out", type=Path, default=None)
ap.add_argument("--tag", default="sensitivity")
a = ap.parse_args()

# float_precision="round_trip" reads every printed digit exactly (the comparison needs it); pandas 3.1 deprecates it
warnings.filterwarnings("ignore", message=".*float_precision.*")
root = Path(os.environ.get("REPRO", "/")).resolve()
show = lambda p: p.resolve().relative_to(root) if p.resolve().is_relative_to(root) else p
KEY = ["flavor", "mass_GeV"]
FLAGS = ["has_sensitivity", "has_sensitivity_N10", "u2_min_open", "u2_max_open", "u2_min_open_N10",
         "u2_max_open_N10"]
txt_new = pd.read_csv(a.new, dtype=str, keep_default_na=False)
txt_ref = pd.read_csv(a.ref, dtype=str, keep_default_na=False)
num_new = pd.read_csv(a.new, float_precision="round_trip")
num_ref = pd.read_csv(a.ref, float_precision="round_trip")
out = []
say = out.append

say(f"rebuilt `{show(a.new)}` ({len(num_new)} rows) vs published `{show(a.ref)}` ({len(num_ref)} rows)")
if list(txt_new.columns) != list(txt_ref.columns):
    say(f"columns differ: {list(txt_new.columns)} vs {list(txt_ref.columns)}")
cols = [c for c in txt_ref.columns if c in txt_new.columns]
for t in (txt_new, txt_ref, num_new, num_ref):
    t["_m"] = pd.to_numeric(t["mass_GeV"]).round(4)
same_order = (len(num_new) == len(num_ref)
              and (num_new.flavor.values == num_ref.flavor.values).all() and (num_new._m.values == num_ref._m.values).all())
say(f"row order (flavour, mass) identical: {same_order}")
j_txt = txt_new.merge(txt_ref, on=["flavor", "_m"], suffixes=("", "_ref"), how="outer", indicator=True)
only = j_txt[j_txt._merge != "both"]
if len(only):
    say(f"rows in one table only: {len(only)}: " + ", ".join(f"{r.flavor} {r._m:g} ({r._merge})" for r in only.itertuples()))
j_txt = j_txt[j_txt._merge == "both"].reset_index(drop=True)
j = num_new.merge(num_ref, on=["flavor", "_m"], suffixes=("", "_ref"))
j = j.set_index(["flavor", "_m"]).loc[list(zip(j_txt.flavor, j_txt._m))].reset_index()
coarse = set(round(m, 4) for m in grids.COARSE)
j["block"] = np.where(j._m.isin(coarse), "coarse", "fine")
j_txt["block"] = j["block"].values
n = len(j)
blocks = ["coarse", "fine"]
nb = {b: int((j.block == b).sum()) for b in blocks}
say(f"matched rows: {n} (coarse {nb['coarse']}, fine {nb['fine']})")
say("")

# 1. identity per column
say("| column | text-identical (all) | coarse | fine | value-identical (all) |")
say("|---|---|---|---|---|")
ident_cols = []
for c in cols:
    if c in KEY:
        continue
    t_same = j_txt[c] == j_txt[f"{c}_ref"]
    x, y = j[c], j[f"{c}_ref"]
    if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
        v_same = (x == y) | (x.isna() & y.isna())
    else:
        v_same = (x.astype(str) == y.astype(str))
    if t_same.all():
        ident_cols.append(c)
    say(f"| {c} | {int(t_same.sum())}/{n} | " + " | ".join(f"{int(t_same[j.block == b].sum())}/{nb[b]}" for b in blocks)
        + f" | {int(v_same.sum())}/{n} |")
say("")
say("bit-identical in every row: " + (", ".join(f"`{c}`" for c in ident_cols) or "none"))
say("")

# 2. relative differences
NUM = ["u2_min", "u2_max", "peak_N", "peak_u2", "u2_min_N10", "u2_max_N10", "peak_u2_N10", "n_events", "n_hits",
       "n_hits_eval"]
NUM = [c for c in NUM if c in cols]
rows = j[["flavor", "_m", "block"]].rename(columns={"_m": "mass_GeV"}).copy()
say("max |new/published - 1| (median in brackets) over rows where both are finite; `k one-sided` = rows where "
    "only one table has a value")
say("")
say("| column | " + " | ".join(f"{fl} {b}" for fl in grids.FLAVORS for b in blocks) + " |")
say("|---|" + "---|" * 6)
for c in NUM:
    x, y = j[c].astype(float), j[f"{c}_ref"].astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where((x == y) | (x.isna() & y.isna()), 0.0, np.abs(x / y - 1))
    both = x.notna() & y.notna() & np.isfinite(r)
    one = x.isna() ^ y.isna()
    rows[f"rel_{c}"] = np.where(both, r, np.nan)
    cells = []
    for fl in grids.FLAVORS:
        for b in blocks:
            s = (j.flavor == fl) & (j.block == b)
            rr = r[s & both]
            cell = f"{rr.max():.2e} ({np.median(rr):.1e})" if len(rr) else "-"
            if (s & one).sum():
                cell += f", {int((s & one).sum())} one-sided"
            cells.append(cell)
    say(f"| {c} | " + " | ".join(cells) + " |")
say("")

# 3. verdicts
flips = []
for c in [f for f in FLAGS if f in cols]:
    x = j[c].astype(str).str.lower() == "true"
    y = j[f"{c}_ref"].astype(str).str.lower() == "true"
    for i in np.flatnonzero(x.values != y.values):
        flips.append((c, j.flavor[i], j._m[i], j.block[i], bool(x.iloc[i]), bool(y.iloc[i]),
                      j.peak_N[i] if "peak_N" in j else np.nan, j.peak_N_ref[i] if "peak_N" in j else np.nan))
if flips:
    say(f"verdict changes: {len(flips)}")
    say("")
    say("| flag | flavour | mass | block | rebuilt | published | peak N rebuilt | peak N published |")
    say("|---|---|---|---|---|---|---|---|")
    for f in flips:
        say(f"| {f[0]} | {f[1]} | {f[2]:g} | {f[3]} | {f[4]} | {f[5]} | {f[6]:.4g} | {f[7]:.4g} |")
else:
    say("verdict changes: none (" + ", ".join(f for f in FLAGS if f in cols) + ")")

text = "\n".join(out)
print(text)
if a.out:
    a.out.mkdir(parents=True, exist_ok=True)
    rows.to_csv(a.out / f"table_{a.tag}_rows.csv", index=False)
    (a.out / f"table_{a.tag}.md").write_text(text + "\n")
