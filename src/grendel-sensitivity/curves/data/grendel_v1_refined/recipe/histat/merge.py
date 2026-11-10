"""Replace the 2.6-3.7 GeV rows of the refined HNL tables with the 4x-statistics scans.

    python merge.py <repo>     # rewrites <repo>/curves/data/grendel_{v1,v2}_refined/hnl/*.csv

v2 tables carry both islands in sensitivity.csv (the N >= 10 file is split from it);
v1 tables keep the v1 layout, so the v1.0 scan output is split first and each island
replaces its own file's rows. Prints how far the window rows moved.
"""
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

H = Path(__file__).resolve().parent
REPO = Path(sys.argv[1])
rd = lambda p: pd.read_csv(p, float_precision="round_trip")
ORDER = {"Ue": 0, "Umu": 1, "Utau": 2}


def split(src: Path, workdir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    workdir.mkdir(exist_ok=True)
    shutil.copyfile(src, workdir / "sensitivity.csv")
    subprocess.run([sys.executable, "-m", "grendel.io.thresholds", str(workdir / "sensitivity.csv"), "--threshold", "10"],
                   check=True, cwd=REPO, env={"PYTHONPATH": str(REPO), "PATH": "/usr/bin:/bin"})
    return rd(workdir / "sensitivity.csv"), rd(workdir / "sensitivity_nsig10.csv")


def replace(table: pd.DataFrame, fresh: pd.DataFrame, label: str) -> pd.DataFrame:
    fresh = fresh[table.columns]
    key = list(zip(fresh.flavor, fresh.mass_GeV.round(6)))
    inwin = [(f, round(m, 6)) in set(key) for f, m in zip(table.flavor, table.mass_GeV)]
    old = table[inwin].merge(fresh, on=["flavor", "mass_GeV"], suffixes=("_old", ""))
    d = (old.u2_min / old.u2_min_old - 1).dropna() * 100
    print(f"  {label}: {len(fresh)} rows replaced; lower edge moved by {d.mean():+.1f} % on average, "
          f"{np.sqrt((d ** 2).mean()):.1f} % rms (expected ~5 % from the old samples' noise)")
    out = pd.concat([table[[not x for x in inwin]], fresh])
    return out.sort_values(["flavor", "mass_GeV"], key=lambda s: s if s.name == "mass_GeV" else s.map(ORDER)).reset_index(drop=True)


for tag in ("v1", "v2"):
    dst = REPO / f"curves/data/grendel_{tag}_refined/hnl"
    s, s10 = rd(dst / "sensitivity.csv"), rd(dst / "sensitivity_nsig10.csv")
    if tag == "v2":
        merged = replace(s, rd(H / "results_v2/sensitivity.csv"), "v2 N>=3+N>=10 table")
        merged.to_csv(dst / "sensitivity.csv", index=False)
        subprocess.run([sys.executable, "-m", "grendel.io.thresholds", str(dst / "sensitivity.csv"), "--threshold", "10"],
                       check=True, cwd=REPO, env={"PYTHONPATH": str(REPO), "PATH": "/usr/bin:/bin"})
    else:
        n3, n10 = split(H / "results_v1/sensitivity.csv", H / "results_v1_split")
        replace(s, n3, "v1 N>=3").to_csv(dst / "sensitivity.csv", index=False)
        replace(s10, n10, "v1 N>=10").to_csv(dst / "sensitivity_nsig10.csv", index=False)
