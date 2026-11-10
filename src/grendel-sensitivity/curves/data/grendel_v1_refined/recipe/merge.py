"""Check the re-scanned 3.2/3.4 GeV rows against the published tables, then write
the refined tables: the published rows plus the new masses.

    python merge.py <v1|v2> <out_dir>
"""
import subprocess
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

W = Path(os.environ.get("REFINE_DIR", Path(__file__).resolve().parent))   # the campaign directory
REPO = Path(__file__).resolve().parents[4]          # this repository
CHECK = (3.2, 3.4)
FIELDS = ("u2_min", "u2_max", "peak_N", "peak_u2")

tag, out = sys.argv[1], Path(sys.argv[2])
published = pd.read_csv(REPO / f"curves/data/grendel_{tag}/hnl/sensitivity.csv", float_precision="round_trip")
new = pd.read_csv(W / f"results_{tag}/sensitivity.csv", float_precision="round_trip")
if list(new.columns) != list(published.columns):
    sys.exit(f"column mismatch:\n  published {list(published.columns)}\n  new       {list(new.columns)}")

worst = 0.0
for _, row in new[new.mass_GeV.round(6).isin(CHECK)].iterrows():
    ref = published[(published.flavor == row.flavor) & np.isclose(published.mass_GeV, row.mass_GeV)].iloc[0]
    for f in FIELDS + tuple(c for c in new.columns if c.endswith("_N10")):
        a, b = float(row[f]), float(ref[f])
        if np.isnan(a) and np.isnan(b):
            continue
        rel = abs(a / b - 1.0) if b else abs(a - b)
        worst = max(worst, rel)
        if rel > 1e-9:
            print(f"  check {tag} {row.flavor} {row.mass_GeV}: {f} new {a:.6e} published {b:.6e} ({rel:.2e})")
print(f"{tag}: re-scanned check rows reproduce the published ones to {worst:.1e}")

added = new[~new.mass_GeV.round(6).isin(CHECK)]
merged = (pd.concat([published, added])
          .sort_values(["flavor", "mass_GeV"], key=lambda s: s if s.name == "mass_GeV" else s.map({"Ue": 0, "Umu": 1, "Utau": 2}))
          .reset_index(drop=True))
dst = out / "hnl"
dst.mkdir(parents=True, exist_ok=True)
merged.to_csv(dst / "sensitivity.csv", index=False)
subprocess.run([sys.executable, "-m", "grendel.io.thresholds", str(dst / "sensitivity.csv"), "--threshold", "10"],
               check=True, cwd=REPO, env={"PYTHONPATH": str(REPO), "PATH": "/usr/bin:/bin"})
pub10 = pd.read_csv(REPO / f"curves/data/grendel_{tag}/hnl/sensitivity_nsig10.csv", float_precision="round_trip")
mine10 = pd.read_csv(dst / "sensitivity_nsig10.csv", float_precision="round_trip")
key = ["flavor", "mass_GeV"]
both = pub10.merge(mine10, on=key, suffixes=("_pub", "_mine"))
cols = [c for c in pub10.columns if c not in key and pd.api.types.is_numeric_dtype(pub10[c])]
same = len(both) == len(pub10) and all(
    np.array_equal(both[f"{c}_pub"].to_numpy(float), both[f"{c}_mine"].to_numpy(float), equal_nan=True) for c in cols)
print(f"{tag}: {len(added)} rows added ({added.mass_GeV.nunique()} masses x {added.flavor.nunique()} flavours); "
      f"published N>=10 rows unchanged by the split: {same}")
