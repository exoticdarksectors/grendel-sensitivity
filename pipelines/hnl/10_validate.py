"""HNL step 10: validate every intermediate and the table against the published campaign outputs.

Writes work/hnl/validation/*.csv and work/hnl/validation/summary.md (also printed). Runs on partial results.

 1. coarse channels   lib/validate_coarse_channels.py: as-built Dmeson/Bc/Kmeson byte for byte; June Bmeson/
                      Bbaryon/induced_tau row by row against the stored combined files; prompt tau and W/Z (MadGraph)
                      weight sums and mean energies.
 2. fine production   summary.json of each scan job against hnl_fine_condor_20261014/results/production_summary.csv:
                      Bmeson, Bc, Bbaryon rows and weight sums exactly (to the last printed digit); W/Z weight sum
                      (= sigma x 1.3) and the MadGraph cross section against results/raw/*/generate_events_run_*.log.
 3. templates         exHad lifetime (deterministic) and, statistically, channel fractions (chi2 of two 20k samples),
                      mean multiplicity, visible fraction against the published MANIFEST records
                      (v2_baseline_20260929/templates/hnl, hnl_fine_condor_20261014/templates_exhad_new).
 4. rows              work/hnl/table against curves/data/grendel_v2_3/hnl: per-mass ratios of the edges and peak N,
                      sensitivity flags, island closure, the Sec. IV numbers; coarse rows also against the lxplus
                      per-job rows.
 5. table             the collected table (step 09) against the published one, row by row (lib/compare_table.py):
                      text-identical columns, max |new/published - 1| per flavour and coarse/fine block, every
                      changed verdict; for sensitivity.csv and sensitivity_nsig10.csv.

Sections 4 and 5 compare with the published table (in this repository) and always run. Sections 1-3 and the lxplus
per-job rows of 4 need the campaign outputs of the optional reference tree $GRENDEL_REF (pipelines/common/ref.py,
reference/README.md); without them they are reported as skipped.
"""
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPRO = Path(os.environ["REPRO"])
sys.path.insert(0, str(REPRO / "pipelines/hnl/lib"))
sys.path.insert(0, str(REPRO / "pipelines/common"))
import grids  # noqa: E402
from ref import PUBLISHED, REF, have_ref  # noqa: E402

W = REPRO / "work/hnl"
V = W / "validation"
V.mkdir(parents=True, exist_ok=True)
CAMP = REF / "grendel/campaigns"
FINE = CAMP / "hnl/hnl_fine_condor_20261014"
PUB = PUBLISHED / "grendel_v2_3/hnl"
rd = lambda p: pd.read_csv(p, float_precision="round_trip")
md = ["# HNL validation (work/hnl/validation)", ""]


def section(title, lines):
    md.extend([f"## {title}", ""] + lines + [""])
    print(f"== {title}\n" + "\n".join(lines), flush=True)


def need(label, *paths):
    """have_ref, plus the line the summary shows when the reference data is missing"""
    if have_ref(label, *paths):
        return None
    rel = [str(p).removeprefix(str(REF) + "/") for p in paths]
    return [f"skipped: needs {', '.join(rel)} of the reference tree (GRENDEL_REF, reference/README.md)"]


# 1. coarse channels
r = subprocess.run([sys.executable, str(REPRO / "pipelines/hnl/lib/validate_coarse_channels.py")],
                   capture_output=True, text=True)
section("1. Coarse-grid channels (step 05, 06)", (r.stdout + r.stderr[-500:]).strip().splitlines())

# 2. fine production
skip2 = need("HNL fine production vs the fine campaign", FINE / "results/production_summary.csv")
ps = rd(FINE / "results/production_summary.csv") if skip2 is None else pd.DataFrame(columns=["flavor", "mass_GeV"])
xs_re = re.compile(r"Cross-section :\s+([0-9.eE+-]+) \+- ([0-9.eE+-]+) pb")
rows = []
for _, p in ps.iterrows():
    fl, m = p.flavor, float(p.mass_GeV)
    s = W / "results/fine" / f"{fl}_{m:g}" / "summary.json"
    lab = grids.label(m)
    d = dict(flavor=fl, mass=m)
    if s.exists():
        j = json.loads(s.read_text())["channels"]
        for ch in ("Bmeson", "Bc", "Bbaryon", "WZ"):
            mine = j.get(ch, {"rows": 0, "w_sum": 0.0})
            ref_n, ref_w = p.get(f"{ch}_rows"), p.get(f"{ch}_w_sum")
            ref_n = 0 if pd.isna(ref_n) else int(ref_n)
            ref_w = 0.0 if pd.isna(ref_w) else float(ref_w)
            d[f"{ch}_rows"], d[f"{ch}_rows_ref"] = mine["rows"], ref_n
            d[f"{ch}_wsum"], d[f"{ch}_wsum_ref"] = mine["w_sum"], ref_w
    lg = W / "fine/wz/llp_4vectors" / fl / "WZ" / f"mN_{lab}.generate_events.log"
    lr = FINE / "results/raw" / f"{fl}_{m:g}" / f"generate_events_run_{lab}.log"
    for key, f in (("xs", lg), ("xs_ref", lr)):
        mm = xs_re.findall(f.read_text()) if f.exists() else []
        if mm:
            d[key], d[key + "_err"] = float(mm[-1][0]), float(mm[-1][1])
    rows.append(d)
f2 = pd.DataFrame(rows)
if skip2 is None:
    f2.to_csv(V / "fine_production.csv", index=False)
lines = []
if "Bmeson_wsum" in f2:
    g = f2.dropna(subset=["Bmeson_wsum"])
    for ch in ("Bmeson", "Bc", "Bbaryon"):
        rows_ok = g[f"{ch}_rows"] == g[f"{ch}_rows_ref"]
        exact = g[f"{ch}_wsum"] == g[f"{ch}_wsum_ref"]
        rel = (g[f"{ch}_wsum"] / g[f"{ch}_wsum_ref"] - 1).abs().replace([np.inf], np.nan).dropna()
        # production_summary.csv prints the Bc sums with 3-4 decimals and the others to 16 digits: compare both ways
        lines.append(f"{ch}: rows identical at {int(rows_ok.sum())}/{len(g)} points; weight sum identical to every "
                     f"printed digit at {int(exact.sum())}, max |new/published - 1| {rel.max() if len(rel) else 0:.1e}")
    w = g[g.WZ_wsum_ref > 0]
    rr = w.WZ_wsum / w.WZ_wsum_ref
    lines.append(f"W/Z weight sum (sigma x 1.3) new/published at {len(w)} points: median {rr.median():.5f}, "
                 f"range {rr.min():.5f}-{rr.max():.5f}")
if "xs" in f2:
    x = f2.dropna(subset=["xs", "xs_ref"])
    # x["xs"], not x.xs: DataFrame.xs is a pandas method
    pull = (x["xs"] - x["xs_ref"]) / np.hypot(x["xs_err"], x["xs_ref_err"])
    lines.append(f"MadGraph cross section at {len(x)} points: ratio median {(x['xs'] / x['xs_ref']).median():.5f}, "
                 f"pull mean {pull.mean():+.2f} rms {np.sqrt((pull ** 2).mean()):.2f}")
section("2. Fine-grid production (scan jobs, step 06)", skip2 or lines or ["no fine points yet"])

# 3. templates
pub = {}
pub_manifests = (CAMP / "v2_baseline_20260929/templates/hnl/MANIFEST.json", FINE / "templates_exhad_new/MANIFEST.json")
skip3 = need("HNL templates vs the published template manifests", *pub_manifests)
for mf in pub_manifests if skip3 is None else ():
    for p in json.loads(mf.read_text())["points"]:
        pub[(p["flavor"], round(p["mass_GeV"], 3))] = p
mine = {}
for mf in sorted((W / "templates_gen").glob("*/MANIFEST.json")):
    for p in json.loads(mf.read_text())["points"]:
        mine[(p["flavor"], round(p["mass_GeV"], 3))] = p
rows = []
N = 20000
for k, a in sorted(mine.items()):
    b = pub.get(k)
    if b is None:
        continue
    fa, fb = a.get("channel_fractions", {}), b.get("channel_fractions", {})
    chans = sorted(set(fa) | set(fb))
    x1 = np.array([fa.get(c, 0.0) for c in chans]); x2 = np.array([fb.get(c, 0.0) for c in chans])
    var = (x1 * (1 - x1) + x2 * (1 - x2)) / N
    use = var > 0
    chi2 = float(((x1 - x2)[use] ** 2 / var[use]).sum())
    rows.append(dict(flavor=k[0], mass=k[1], ctau=a["ctau_m_u2eq1"], ctau_ref=b["ctau_m_u2eq1"],
                     chi2_channels=chi2, ndf=int(use.sum()) - 1,
                     mult=a["mean_multiplicity"], mult_ref=b["mean_multiplicity"],
                     vis=a["vis_frac_ge2_charged"], vis_ref=b["vis_frac_ge2_charged"]))
t3 = pd.DataFrame(rows)
lines = []
if len(t3):
    t3.to_csv(V / "templates.csv", index=False)
    from scipy.stats import chi2 as C2
    pv = C2.sf(t3.chi2_channels, t3.ndf.clip(lower=1))
    lines += [f"{len(t3)} templates compared",
              f"ctau(U^2=1) identical at {int((t3.ctau == t3.ctau_ref).sum())}/{len(t3)}; max rel diff "
              f"{(t3.ctau / t3.ctau_ref - 1).abs().max():.1e}",
              f"channel fractions: chi2/ndf summed {t3.chi2_channels.sum():.0f}/{t3.ndf.sum()}; p-values: "
              f"{int((pv < 0.01).sum())} below 0.01 (expected {0.01 * len(t3):.1f})",
              f"mean multiplicity ratio {(t3.mult / t3.mult_ref).median():.4f} "
              f"({(t3.mult / t3.mult_ref).min():.4f}-{(t3.mult / t3.mult_ref).max():.4f}); visible fraction diff "
              f"{(t3.vis - t3.vis_ref).mean():+.4f} (rms {np.sqrt(((t3.vis - t3.vis_ref) ** 2).mean()):.4f}, "
              f"binomial sigma ~{np.sqrt(2 * 0.8 * 0.2 / N):.4f})"]
section("3. Templates (step 07)", skip3 or lines or ["no templates yet"])

# 4. rows
lines = []
mine_rows = []
for g, ms in (("coarse", grids.COARSE), ("fine", grids.FINE)):
    for fl in grids.FLAVORS:
        for m in ms:
            p = W / "results" / g / f"{fl}_{m:g}" / "sensitivity.csv"
            if p.exists():
                mine_rows.append(rd(p).assign(grid=g))
if mine_rows:
    t = pd.concat(mine_rows, ignore_index=True)
    ref = rd(PUB / "sensitivity.csv")
    j = t.merge(ref, on=["flavor", "mass_GeV"], suffixes=("", "_ref"))
    for c in ("u2_min", "u2_max", "peak_N", "u2_min_N10", "u2_max_N10"):
        j[f"r_{c}"] = j[c] / j[f"{c}_ref"]
    j.to_csv(V / "rows.csv", index=False)
    lines.append(f"{len(j)} of 723 rows; sensitivity flag (N>=3) agrees at "
                 f"{int((j.has_sensitivity.astype(bool) == j.has_sensitivity_ref.astype(bool)).sum())}, "
                 f"(N>=10) at {int((j.has_sensitivity_N10.astype(bool) == j.has_sensitivity_N10_ref.astype(bool)).sum())}")
    lines.append("| flavour | mass range | rows | u2_min new/pub median (min-max) | u2_max | peak N |")
    lines.append("|---|---|---|---|---|---|")
    for fl in grids.FLAVORS:
        for lo, hi in ((0, 0.5), (0.5, 2.0), (2.0, 3.4), (3.4, 4.2), (4.2, 11)):
            s = j[(j.flavor == fl) & (j.mass_GeV >= lo) & (j.mass_GeV < hi)]
            if not len(s):
                continue
            f = lambda c: (lambda x: f"{x.median():.3f} ({x.min():.3f}-{x.max():.3f})" if len(x) else "-")(
                s[c].replace([np.inf, -np.inf], np.nan).dropna())
            lines.append(f"| {fl} | {lo}-{hi} | {len(s)} | {f('r_u2_min')} | {f('r_u2_max')} | {f('r_peak_N')} |")
    LX = CAMP / "v2_3_20261016/lxplus/out"
    if (j.grid == "coarse").any() and (skip4 := need("HNL coarse rows vs the lxplus job rows", LX)) is not None:
        lines += skip4
    elif (j.grid == "coarse").any():
        lx = []
        for _, rrow in j[j.grid == "coarse"].iterrows():
            p = LX / f"hnl_{rrow.flavor}_{rrow.mass_GeV:g}" / "sensitivity.csv"
            if p.exists():
                lx.append(rd(p))
        if lx:
            lxt = pd.concat(lx)
            jj = t.merge(lxt, on=["flavor", "mass_GeV"], suffixes=("", "_lx"))
            same = sum(np.array_equal(a, b) for a, b in zip(jj.u2_min.values, jj.u2_min_lx.values))
            lines.append(f"coarse rows vs the lxplus job rows: {len(jj)} compared, u2_min bit-identical at {same}")
    # closure and Sec. IV numbers (full table only)
    for fl in grids.FLAVORS:
        a = j[(j.flavor == fl)]
        sens = a[a.has_sensitivity.astype(bool)]
        sref = a[a.has_sensitivity_ref.astype(bool)]
        if len(sens):
            best = sens.loc[sens.u2_min.idxmin()]
            bref = sref.loc[sref.u2_min_ref.idxmin()]
            lines.append(f"{fl}: last sensitive mass {sens.mass_GeV.max():g} GeV (published {sref.mass_GeV.max():g}); "
                         f"best |U|^2 {best.u2_min:.3g} at {best.mass_GeV:g} GeV (published {bref.u2_min_ref:.3g} at "
                         f"{bref.mass_GeV:g})")
    for fl, m, quoted in (("Umu", 1.6, "1.3e-8"), ("Ue", 1.6, "1.3e-8"), ("Ue", 0.395, "1.1e-8 near 0.4")):
        s = j[(j.flavor == fl) & (np.isclose(j.mass_GeV, m))]
        if len(s):
            lines.append(f"Sec. IV: {fl} at {m:g} GeV u2_min {s.u2_min.iloc[0]:.3g} (published row "
                         f"{s.u2_min_ref.iloc[0]:.3g}; text {quoted})")
    if len(j[(j.flavor == "Ue") & (j.mass_GeV < 0.5)]):
        s = j[(j.flavor == "Ue") & (j.mass_GeV < 0.5)]
        lines.append(f"Sec. IV: Ue min below 0.5 GeV {s.u2_min.min():.3g} at {s.loc[s.u2_min.idxmin()].mass_GeV:g} "
                     f"(published {s.u2_min_ref.min():.3g}; text 1.1e-8 near 0.4 GeV)")
    if len(j[j.flavor == "Utau"]):
        s = j[(j.flavor == "Utau") & j.has_sensitivity.astype(bool)]
        if len(s):
            lines.append(f"Sec. IV: Utau best {s.u2_min.min():.3g} at {s.loc[s.u2_min.idxmin()].mass_GeV:g} GeV "
                         f"(text ~1.5e-7)")
section("4. Rows (step 08, 09)", lines or ["no rows yet"])

# 5. the table
T = W / "table/hnl"
lines = []
for name in ("sensitivity", "sensitivity_nsig10"):
    if (T / f"{name}.csv").exists():
        r = subprocess.run([sys.executable, str(REPRO / "pipelines/hnl/lib/compare_table.py"), str(T / f"{name}.csv"),
                            str(PUB / f"{name}.csv"), "--out", str(V), "--tag", name], capture_output=True, text=True)
        lines += [f"### {name}.csv", ""] + (r.stdout + r.stderr[-500:]).strip().splitlines() + [""]
section("5. Table (step 09) against the published table", lines or ["no table yet (step 09)"])
(V / "summary.md").write_text("\n".join(md) + "\n")
