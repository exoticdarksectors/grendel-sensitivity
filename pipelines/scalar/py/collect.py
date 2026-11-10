"""Build the BC4 and BC5 tables from the one-mass scan folders, as scalar_eventcalc_20261019/collect.py does.

    python collect.py WORK RUN RESULTS POINTS_BC4 POINTS_BC5 SRC

WORK/<bench>/scans/RUN/out/<mass>/sensitivity.csv -> RESULTS/<bench>/sensitivity.csv (one row per mass, sorted by
mass, pandas to_csv as the original), sensitivity_nsig10.csv (grendel.io.thresholds --threshold 10), and for BC5
run_metadata.json (merged from the per-job files with the published recipe's fields: workers "1 per job", code
c621279 = the scan code of 8802f50, geometry 2026-10-16.1; config paths relative to the work folder).
Missing masses are listed and leave the table incomplete (exit status 1).
"""
import json
import subprocess
import sys
from pathlib import Path

import warnings

import pandas as pd

warnings.filterwarnings("ignore", message=".*float_precision.*")

work, run, results = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
points = {"bc4": sys.argv[4].split(), "bc5": sys.argv[5].split()}
rd = lambda p: pd.read_csv(p, float_precision="round_trip")
status = 0
for bench in ("bc4", "bc5"):
    out = work / bench / "scans" / run / "out"
    rows, missing = [], []
    for m in points[bench]:
        p = out / m / "sensitivity.csv"
        (rows.append(rd(p)) if p.exists() else missing.append(m))
    if not rows:
        print(f"{bench}: no rows yet")
        status = 1
        continue
    t = pd.concat(rows, ignore_index=True).sort_values("mass_GeV").reset_index(drop=True)
    d = results / bench
    d.mkdir(parents=True, exist_ok=True)
    t.to_csv(d / "sensitivity.csv", index=False)
    subprocess.run([sys.executable, "-m", "grendel.io.thresholds", str(d / "sensitivity.csv"), "--threshold", "10"],
                   check=True, capture_output=True)
    print(f"{bench}: {len(t)}/{len(points[bench])} masses -> {d}" + (f"; MISSING {' '.join(missing)}" if missing else ""))
    if missing:
        status = 1
    if bench == "bc5":
        meta = [json.loads((out / m / "run_metadata.json").read_text()) for m in points[bench]
                if (out / m / "run_metadata.json").exists()]
        m0 = {k: v for k, v in meta[0].items() if k not in ("timestamp", "masses_GeV", "n_points_requested",
                                                            "n_points_processed", "n_points_skipped", "skipped")}
        m0.update(masses_GeV=sorted(t.mass_GeV.tolist()), n_points_requested=len(points[bench]),
                  n_points_processed=len(t), workers="1 per job (pipelines/scalar/scan_one.sh)", code="c621279",
                  geometry="grendel-geometry 2026-10-16.1")
        for k in [k for k in m0 if k.startswith("n_sensitive")]:
            m0.pop(k)
        m0["n_sensitive"] = int(t.has_sensitivity.astype(bool).sum())
        m0["inputs"] = {k: v for m in meta for k, v in m.get("inputs", {}).items()}
        m0["total_time_s"] = round(sum(m.get("total_time_s", 0) for m in meta), 1)
        for k in ("vectors", "templates", "templates_dir"):
            if k in m0.get("config", {}).get("model", {}):
                m0["config"]["model"][k] = str(Path(m0["config"]["model"][k]).relative_to(work.parent.parent)) \
                    if str(m0["config"]["model"][k]).startswith(str(work.parent.parent)) else m0["config"]["model"][k]
        (d / "run_metadata.json").write_text(json.dumps(m0, indent=2) + "\n")
sys.exit(status)
