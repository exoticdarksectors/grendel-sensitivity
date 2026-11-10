"""v1 variant of merge.py. The published v1 tables are in the v1 layout: no _N10
columns in sensitivity.csv, the N >= 10 island in its own sensitivity_nsig10.csv
with the same 14 columns. The v1.0 scan writes both islands in one table, so the
check compares each island against its own published file, and the new rows are
split with grendel.io.thresholds before they are appended.

    python merge_v1.py <out_dir>
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
out = Path(sys.argv[1])
rd = lambda p: pd.read_csv(p, float_precision="round_trip")

pub = rd(REPO / "curves/data/grendel_v1/hnl/sensitivity.csv")
pub10 = rd(REPO / "curves/data/grendel_v1/hnl/sensitivity_nsig10.csv")
scratch = W / "results_v1_split"
scratch.mkdir(exist_ok=True)
(scratch / "sensitivity.csv").write_bytes((W / "results_v1/sensitivity.csv").read_bytes())
subprocess.run([sys.executable, "-m", "grendel.io.thresholds", str(scratch / "sensitivity.csv"), "--threshold", "10"],
               check=True, cwd=REPO, env={"PYTHONPATH": str(REPO), "PATH": "/usr/bin:/bin"})
new = rd(scratch / "sensitivity.csv")[pub.columns]
new10 = rd(scratch / "sensitivity_nsig10.csv")[pub10.columns]

worst = 0.0
for table, ref_table, label in ((new, pub, "N>=3"), (new10, pub10, "N>=10")):
    for _, row in table[table.mass_GeV.round(6).isin(CHECK)].iterrows():
        ref = ref_table[(ref_table.flavor == row.flavor) & np.isclose(ref_table.mass_GeV, row.mass_GeV)].iloc[0]
        for f in FIELDS:
            a, b = float(row[f]), float(ref[f])
            if np.isnan(a) and np.isnan(b):
                continue
            rel = abs(a / b - 1.0)
            worst = max(worst, rel)
            if rel > 1e-9:
                print(f"  check v1 {label} {row.flavor} {row.mass_GeV}: {f} new {a:.6e} published {b:.6e} ({rel:.2e})")
print(f"v1: re-scanned check rows reproduce the published ones to {worst:.1e}")

order = {"Ue": 0, "Umu": 1, "Utau": 2}
dst = out / "hnl"
dst.mkdir(parents=True, exist_ok=True)
for name, published, fresh in (("sensitivity.csv", pub, new), ("sensitivity_nsig10.csv", pub10, new10)):
    added = fresh[~fresh.mass_GeV.round(6).isin(CHECK)]
    merged = (pd.concat([published, added])
              .sort_values(["flavor", "mass_GeV"], key=lambda s: s if s.name == "mass_GeV" else s.map(order))
              .reset_index(drop=True))
    merged.to_csv(dst / name, index=False)
    print(f"v1 {name}: {len(added)} rows added")
