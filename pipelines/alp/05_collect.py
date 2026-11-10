"""Step 05 -- build the BC10 tables from the per-mass scans and compare them with the published v2.4 table.

As bc10_v24_20261019/collect.py: one row per job, sorted by mass, written with pandas (shortest round-trip floats),
then `python -m grendel.io.thresholds sensitivity.csv --threshold 10` for sensitivity_nsig10.csv.

    python 05_collect.py VARIANT         VARIANT = main | control (work/alp/scan_<VARIANT>/bc10_<M>/)

Outputs: work/alp/tables_<VARIANT>/bc10/sensitivity{,_nsig10}.csv and work/alp/validation/rows_<VARIANT>.txt with
  * the row-by-row comparison against curves/data/grendel_v2_4/bc10 (bit-identical fields, n_hits, largest relative
    difference of peak_N and of the edges at N >= 3 and N >= 10, has_sensitivity flips);
  * the published checks: sensitive masses (N >= 3), reach, deepest point (1/f and g_Y = v_h/f, v_h = 246 GeV) and
    the lower edge below 1 GeV, as quoted in the paper (g_Y ~ 2.5e-6 at 1.2-1.5 GeV, (0.6-2.0)e-5 below 1 GeV).
Exits non-zero only when rows are missing.
"""
import hashlib
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=".*float_precision.*")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from ref import PUBLISHED, REF, have_ref  # noqa: E402  (pipelines/common/ref.py)

REPRO = Path(__file__).resolve().parents[2]
W = REPRO / "work/alp"
PUB = PUBLISHED / "grendel_v2_4/bc10"          # the published table (in this repository)
STORED_TEMPLATES = REF / "grendel/campaigns/bc10/bc10_v24_20261019/templates_exhad"   # control variant ($GRENDEL_REF)
POINTS = Path(__file__).resolve().parent / "points_all.txt"
VH = 246.0
EDGES = ["peak_N", "invf_min", "invf_max", "invf_min_N10", "invf_max_N10", "peak_invf"]


def rd(path):
    return pd.read_csv(path, float_precision="round_trip")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def summary(t, tag):
    s = t[t.has_sensitivity.astype(bool)]
    best = s.loc[s.invf_min.idxmin()]
    low = s[s.mass_GeV < 1.0]
    deep = s[s.mass_GeV.between(1.2, 1.5)]
    return [f"{tag}: sensitive (N>=3) {len(s)}, (N>=10) {int(t.has_sensitivity_N10.astype(bool).sum())}; "
            f"reach (N>=3) {s.mass_GeV.max():g} GeV; deepest 1/f = {best.invf_min:.4g} GeV^-1 "
            f"(g_Y = {VH * best.invf_min:.3g}) at {best.mass_GeV:g} GeV; g_Y lower edge at 1.2-1.5 GeV "
            f"{VH * deep.invf_min.min():.3g}-{VH * deep.invf_min.max():.3g}, below 1 GeV "
            f"{VH * low.invf_min.min():.3g}-{VH * low.invf_min.max():.3g}"]


def main(variant):
    scans = W / f"scan_{variant}"
    masses = POINTS.read_text().split()
    rows, missing = [], []
    for m in masses:
        p = scans / f"bc10_{m}" / "sensitivity.csv"
        (rows.append(rd(p)) if p.exists() else missing.append(m))
    if not rows:
        if variant == "control" and not have_ref("BC10 control table (needs the stored templates)", STORED_TEMPLATES):
            return 0
        print("no rows yet")
        return 1
    t = pd.concat(rows, ignore_index=True).sort_values("mass_GeV").reset_index(drop=True)
    d = W / f"tables_{variant}/bc10"
    d.mkdir(parents=True, exist_ok=True)
    t.to_csv(d / "sensitivity.csv", index=False)
    subprocess.run([sys.executable, "-m", "grendel.io.thresholds", str(d / "sensitivity.csv"), "--threshold", "10"],
                   check=True, capture_output=True)

    pub = rd(PUB / "sensitivity.csv")
    rep = [f"# BC10 {variant} scan against curves/data/grendel_v2_4/bc10", "",
           f"rows {len(t)} of {len(masses)}" + (f"; MISSING {len(missing)}: {' '.join(missing)}" if missing else ""),
           f"sha256 sensitivity.csv        {sha(d / 'sensitivity.csv')}  (published {sha(PUB / 'sensitivity.csv')})",
           f"sha256 sensitivity_nsig10.csv {sha(d / 'sensitivity_nsig10.csv')}  "
           f"(published {sha(PUB / 'sensitivity_nsig10.csv')})"]
    m = t.merge(pub, on="mass_GeV", suffixes=("", "_pub"))
    cols = [c for c in t.columns if c != "mass_GeV"]

    def same(a, b):
        return (a == b) | (pd.isna(a) & pd.isna(b))

    identical = np.all([same(m[c], m[c + "_pub"]) for c in cols], axis=0)
    rep.append(f"rows bit-identical to the published row in every field: {int(identical.sum())} / {len(m)}")
    rep.append(f"n_hits identical: {int((m.n_hits == m.n_hits_pub).sum())} / {len(m)}; "
               f"has_sensitivity flips (N>=3): {int((~same(m.has_sensitivity, m.has_sensitivity_pub)).sum())}, "
               f"(N>=10): {int((~same(m.has_sensitivity_N10, m.has_sensitivity_N10_pub)).sum())}"
               " (empty == empty: the N>=10 flag is empty at the masses inside the vetoed resonance windows)")
    for c in EDGES:
        r = (m[c] / m[c + "_pub"] - 1).abs().replace([np.inf], np.nan).dropna()
        if len(r):
            worst = m.loc[r.idxmax(), "mass_GeV"]
            rep.append(f"  {c:14s}: max |rel diff| {r.max():.2e} (at {worst:g} GeV), median {r.median():.1e}, "
                       f"> 1e-3 at {int((r > 1e-3).sum())} masses, > 1e-2 at {int((r > 1e-2).sum())}")
    n10 = rd(d / "sensitivity_nsig10.csv").merge(rd(PUB / "sensitivity_nsig10.csv"), on="mass_GeV", suffixes=("", "_pub"))
    c10 = [c for c in n10.columns if not c.endswith("_pub") and c != "mass_GeV" and c + "_pub" in n10.columns]
    rel = [(0.0, "", float("nan"))]
    for c in c10:
        if pd.api.types.is_numeric_dtype(n10[c]) and not pd.api.types.is_bool_dtype(n10[c]):
            x = (n10[c] / n10[c + "_pub"] - 1).abs().replace([np.inf], np.nan).dropna()
            if len(x):
                rel.append((x.max(), c, n10.loc[x.idxmax(), "mass_GeV"]))
    big = max(rel, key=lambda r: r[0])
    rep.append(f"sensitivity_nsig10.csv: largest |rel diff| over all numeric fields {big[0]:.1e} ({big[1]} at {big[2]:g} GeV)")
    rep += [""] + summary(t, "this run") + summary(pub, "published")
    rep += ["", "per-mass ratios (this run / published), masses where any edge differs by more than 1e-6:"]
    for _, r in m.iterrows():
        diffs = {c: r[c] / r[c + "_pub"] for c in EDGES
                 if np.isfinite(r[c]) and np.isfinite(r[c + "_pub"]) and r[c + "_pub"] != 0
                 and abs(r[c] / r[c + "_pub"] - 1) > 1e-6}
        if diffs or r.has_sensitivity != r.has_sensitivity_pub:
            rep.append(f"  {r.mass_GeV:g}: n_hits {r.n_hits}/{r.n_hits_pub} "
                       + " ".join(f"{k} {v:.6f}" for k, v in diffs.items())
                       + ("" if r.has_sensitivity == r.has_sensitivity_pub else
                          f" has_sensitivity {r.has_sensitivity}/{r.has_sensitivity_pub}"))
    out = W / f"validation/rows_{variant}.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(rep) + "\n")
    print("\n".join(rep[:16 + len(EDGES)]))
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "main"))
