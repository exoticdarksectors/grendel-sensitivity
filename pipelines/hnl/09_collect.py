"""HNL step 09: the table from the per-point rows, as v2_3_20261016/collect.py built it.

One row per job, the 104 x 3 coarse (lxplus) points then the 137 x 3 fine points, sorted by flavour (Ue, Umu, Utau)
then mass; read with float_precision="round_trip" and written with pandas' to_csv(index=False); the N >= 10 island
split off with grendel.io.thresholds --threshold 10. Writes work/hnl/table/hnl/sensitivity{,_nsig10}.csv (a
results tree for curves.plot). Exits 1 while points are missing (it lists them).
"""
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPRO = Path(os.environ["REPRO"])
sys.path.insert(0, str(REPRO / "pipelines/hnl/lib"))
import grids  # noqa: E402

R = REPRO / "work/hnl/results"
ORDER = {"Ue": 0, "Umu": 1, "Utau": 2}
rd = lambda p: pd.read_csv(p, float_precision="round_trip")
rows, missing = [], []
for g, ms in (("coarse", grids.COARSE), ("fine", grids.FINE)):
    for fl in grids.FLAVORS:
        for m in ms:
            p = R / g / f"{fl}_{m:g}" / "sensitivity.csv"
            (rows.append(rd(p)) if p.exists() else missing.append(f"{g} {fl} {m:g}"))
if missing:
    print(f"{len(missing)} of 723 points missing, e.g. {missing[:10]}")
    sys.exit(1)
t = pd.concat(rows, ignore_index=True)
t = t.sort_values(["flavor", "mass_GeV"], key=lambda s: s if s.name == "mass_GeV" else s.map(ORDER)).reset_index(drop=True)
d = REPRO / "work/hnl/table/hnl"
d.mkdir(parents=True, exist_ok=True)
t.to_csv(d / "sensitivity.csv", index=False)
subprocess.run([sys.executable, "-m", "grendel.io.thresholds", str(d / "sensitivity.csv"), "--threshold", "10"],
               check=True, capture_output=True)
print(f"table: {len(t)} rows -> {d}")
