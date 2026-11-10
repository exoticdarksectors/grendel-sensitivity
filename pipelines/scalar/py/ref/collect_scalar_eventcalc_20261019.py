"""Build the BC4 and BC5 tables with the EventCalc scalar model from the condor outputs; compare with v2.3 (Winkler).

    rsync -a <lxplus>:<campaign work dir>/grendel_scalar_ec/out/ lxplus/out/
    PYTHONPATH=../../bc10/bc10_v24_20261019/repo_private python collect.py

tables/<bc>/sensitivity.csv (+ _nsig10 from grendel.io.thresholds): one row per job, sorted by mass.
comparison/summary.md: sensitive masses, reach, deepest point and edges new/v2.3 per mass range, next to the lifetime
ratio EventCalc/Winkler, so the edge shifts can be read against the lifetime change; then the Sec. IV numbers
(v2_3_20261016/edge_ratios.py for BC4, bc5_v2_1_20261003/results/bc5_hybrid/paper_numbers.py for BC5).
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

C = Path(__file__).resolve().parent
CAMP = C.parents[1]
REPO = CAMP / "bc10/bc10_v24_20261019/repo_private"
OUT = C / "lxplus/out"
RANGES = {"bc4": [(0.2, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 3.5), (3.5, 5.0)],
          "bc5": [(0.2, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 3.5), (3.5, 5.0), (5.0, 10.0), (10.0, 63.0)]}
rd = lambda p: pd.read_csv(p, float_precision="round_trip")
EC = np.array(__import__("json").load(__import__("gzip").open(
    REPO / "grendel/models/scalar/data/eventcalc_1809/ctau.json.gz")), float)


def ctau_ec(m):
    """c tau (sin^2 theta = 1) of the EventCalc table, the lifetime of the exHad scalar-1809 templates."""
    return float(np.exp(np.interp(np.log(m), np.log(EC[:, 0]), np.log(EC[:, 1]))))


def masses(bench):
    return [l.split()[2] for l in (C / f"lxplus/points_{bench}.txt").read_text().splitlines() if l.strip()]


def build(bench):
    rows, missing = [], []
    for m in masses(bench):
        p = OUT / f"{bench}_-_{m}" / "sensitivity.csv"
        (rows.append(rd(p)) if p.exists() else missing.append(m))
    t = pd.concat(rows, ignore_index=True).sort_values("mass_GeV").reset_index(drop=True)
    d = C / "tables" / bench
    d.mkdir(parents=True, exist_ok=True)
    t.to_csv(d / "sensitivity.csv", index=False)
    subprocess.run([sys.executable, "-m", "grendel.io.thresholds", str(d / "sensitivity.csv"), "--threshold", "10"],
                   check=True, capture_output=True)
    return t, missing


def bc5_metadata(t):
    """One run_metadata.json for the BC5 table from the per-job files (as v2_3_20261016/collect.py)."""
    import json
    meta = [json.loads(p.read_text()) for p in sorted(OUT.glob("bc5_-_*/run_metadata.json"))]
    m0 = {k: v for k, v in meta[0].items() if k not in ("timestamp", "masses_GeV", "n_points_requested",
                                                        "n_points_processed", "n_points_skipped", "skipped")}
    m0.update(masses_GeV=sorted(t.mass_GeV.tolist()), n_points_requested=len(masses("bc5")),
              n_points_processed=len(t), workers="1 per HTCondor job", code="c621279",
              geometry="grendel-geometry 2026-10-16.1")
    for k in [k for k in m0 if k.startswith("n_sensitive")]:
        m0.pop(k)
    m0["n_sensitive"] = int(t.has_sensitivity.astype(bool).sum())
    m0["inputs"] = {k: v for m in meta for k, v in m.get("inputs", {}).items()}
    m0["total_time_s"] = round(sum(m.get("total_time_s", 0) for m in meta), 1)
    for k in ("vectors", "templates", "templates_dir"):          # lxplus paths -> the inputs they hold
        if k in m0.get("config", {}).get("model", {}):
            m0["config"]["model"][k] = m0["config"]["model"][k].split("/grendel_scalar_ec/")[-1]
    (C / "tables/bc5/run_metadata.json").write_text(json.dumps(m0, indent=2) + "\n")


def compare(bench, t, missing):
    from grendel.models.scalar import model
    old = rd(REPO / f"curves/data/grendel_v2_3/{bench}/sensitivity.csv")
    sens = lambda x: x[x.has_sensitivity.astype(bool)]
    s, so = sens(t), sens(old)
    rep = [f"## {bench.upper()}: EventCalc scalar model / v2.3 (Winkler hybrid)", "",
           f"rows {len(t)}" + (f", MISSING {len(missing)}: {' '.join(missing)}" if missing else "") +
           f"; sensitive (N>=3) {len(s)} / {len(so)}; N>=10 {int(t.has_sensitivity_N10.sum())} / "
           f"{int(old.has_sensitivity_N10.sum())}",
           f"reach (last sensitive mass, N>=3): {s.mass_GeV.max():g} / {so.mass_GeV.max():g} GeV",
           f"deepest lower edge: {s.u2_min.min():.3g} at {s.loc[s.u2_min.idxmin(), 'mass_GeV']:g} GeV / "
           f"{so.u2_min.min():.3g} at {so.loc[so.u2_min.idxmin(), 'mass_GeV']:g} GeV", "",
           "edges new/v2.3 per mass range, median (min-max), and c tau EventCalc/Winkler there:"]
    m = t.merge(old, on="mass_GeV", suffixes=("", "_old"))
    m = m[m.has_sensitivity.astype(bool) & m.has_sensitivity_old.astype(bool)]
    for lo, hi in RANGES[bench]:
        x = m[m.mass_GeV.between(lo, hi - 1e-9)]
        if not len(x):
            continue
        r1, r2 = x.u2_min / x.u2_min_old, (x.u2_max / x.u2_max_old).dropna()
        lt = np.array([ctau_ec(v) / model.ctau_sin2theta1(v, "winkler") for v in x.mass_GeV])
        rep.append(f"  {lo}-{hi} GeV ({len(x)}): lower {r1.median():.3f} ({r1.min():.3f}-{r1.max():.3f}); "
                   f"upper {r2.median():.3f} ({r2.min():.3f}-{r2.max():.3f})"
                   + f"; ctau EventCalc/Winkler {np.median(lt):.2f} ({lt.min():.2f}-{lt.max():.2f})")
    return rep


def main():
    rep = ["# BC4 and BC5 with the EventCalc scalar model (exHad scalar-1809 lifetime, BRs and final states)", ""]
    for bench in ("bc4", "bc5"):
        t, missing = build(bench)
        if bench == "bc5":
            bc5_metadata(t)
        rep += compare(bench, t, missing) + [""]
    env = dict(__import__("os").environ, PYTHONPATH=str(REPO))
    rep += ["## Sec. IV numbers", "", "### BC4 (v2_3_20261016/edge_ratios.py)", "```",
            subprocess.run([sys.executable, str(CAMP / "v2_3_20261016/edge_ratios.py"), "bc4",
                            str(C / "tables/bc4/sensitivity.csv")], capture_output=True, text=True, env=env).stdout,
            "```", "", "### BC5 (bc5_hybrid/paper_numbers.py)", "```",
            subprocess.run([sys.executable, str(CAMP / "bc5/bc5_v2_1_20261003/results/bc5_hybrid/paper_numbers.py"),
                            str(C / "tables/bc5/sensitivity.csv"), str(REPO)], capture_output=True, text=True,
                           env=env).stdout, "```"]
    (C / "comparison").mkdir(exist_ok=True)
    (C / "comparison/summary.md").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
