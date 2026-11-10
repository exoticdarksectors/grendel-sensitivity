"""Scan one (flavour, mass) point of the v2.3 HNL table, as one lxplus job did.

    python scan_point.py coarse|fine <flavor> <mass>

coarse  the combined file of the as-built 1x production, rebuilt from the kept channel files (step 05: Dmeson, Bc,
        Kmeson as-built; Bmeson, Bbaryon, tau, induced_tau June v1; step 06: W/Z 100k), in the order of
        combine_channels.CHANNELS (the order of the stored asbuilt_20261014 combined files);
        scan: v2_3_20261016/lxplus/job_v23.sh.
fine    hnl_fine_condor_20261014/condor/job.sh: B, B_c, Lambda_b pools of 1 000 000 (run_all --channels bottom bc
        baryon, seed 7000/7100/7200, this mass alone) regenerated here (seconds; deterministic), the W/Z 1M sample
        of step 06, the closed channels as empty files, combine_for_point(strict) over Bmeson, Dmeson, Bc, Bbaryon,
        induced_tau, Kmeson, WZ; scan: v2_3_20261016/hnl_fine/condor/job_scan.sh.
Both: grendel.models.hnl.scan --flavor F --mass M --workers 1 --thresholds 3 10 (100 decay samples), templates of
step 07, a per-job geometry cache. Writes work/hnl/results/<grid>/<flavor>_<mass>/ (sensitivity.csv, scan.log,
summary.json with rows, weight sums and SHA-256 per channel and of the combined file); the combined file and the
job's channel files are deleted afterwards. Resume-safe: a point with sensitivity.csv is skipped.
A fine job (four 1M pools, 3.2 GB peak) first passes lib/memgate.py: at most HNL_FINE_SCAN_MAX (2) at once, and only
when the memory is free.
"""
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

grid, flavor, mass = sys.argv[1], sys.argv[2], float(sys.argv[3])
REPRO = Path(os.environ["REPRO"])
sys.path.insert(0, str(REPRO / "pipelines/hnl/lib"))
import grids  # noqa: E402

W = REPRO / "work/hnl"
lab = grids.label(mass)
mtxt = f"{mass:g}"
out = W / "results" / grid / f"{flavor}_{mtxt}"
if (out / "sensitivity.csv").exists():
    print(f"done: {out}")
    sys.exit(0)
if grid == "fine":
    import memgate  # noqa: E402
    _token = memgate.acquire(log=lambda s: print(s, flush=True))      # held until this process ends
job = W / "jobs" / f"{grid}_{flavor}_{mtxt}"
shutil.rmtree(job, ignore_errors=True)
vec = job / "v"
(vec / flavor).mkdir(parents=True)
env = dict(os.environ, GRENDEL_HNL_VECTORS_DIR=str(vec), GRENDEL_HNL_CACHE_DIR=str(job / "cache"),
           GRENDEL_WORK_DIR=str(job / "work"), MPLCONFIGDIR=str(job / "mpl"), PYTHONDONTWRITEBYTECODE="1")
PROD = str(W / "venv_prod/bin/python")
t0 = time.time()
log = lambda s: print(f"[{time.time() - t0:7.1f}s] {s}", flush=True)


def place(src: Path, ch: str):
    dst = vec / flavor / ch / f"mN_{lab}.csv"
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.suffix == ".gz":
        with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
            shutil.copyfileobj(fi, fo)
    else:
        os.symlink(src, dst)


if grid == "coarse":
    C = W / "coarse"
    for prod, chans in (("asbuilt", ["Dmeson", "Bc", "Kmeson"]), ("v1", ["Bmeson", "Bbaryon", "tau", "induced_tau"]),
                        ("wz", ["WZ"])):
        for ch in chans:
            src = C / prod / "llp_4vectors" / flavor / ch / f"mN_{lab}.csv"
            if not src.exists():
                sys.exit(f"missing input {src}")
            place(src, ch)
    channels = None                                   # combine_channels.CHANNELS
elif grid == "fine":
    seed = grids.SEEDS_FINE[flavor]
    r = subprocess.run([PROD, "-m", "grendel.models.hnl.production.run_all", "--flavor", flavor, "--masses", mtxt,
                        "--channels", "bottom", "bc", "baryon", "--no-wz", "--no-prompt-tau", "--skip-combine",
                        "--n-pool", "1000000", "--workers", "1", "--seed", str(seed)],
                       env=env, cwd=job, stdout=open(job / "prod_mesons.log", "w"), stderr=subprocess.STDOUT)
    if r.returncode:
        sys.exit(f"meson production failed ({job / 'prod_mesons.log'})")
    log("mesons produced")
    src = W / "fine/wz/llp_4vectors" / flavor / "WZ" / f"mN_{lab}.csv.gz"
    if not src.exists():
        sys.exit(f"missing input {src}")
    place(src, "WZ")
    for ch in ("Dmeson", "induced_tau", "Kmeson", "tau", "Bbaryon"):      # closed above 2 GeV, as published
        p = vec / flavor / ch / f"mN_{lab}.csv"
        if not p.exists():
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text("")
    channels = ["Bmeson", "Dmeson", "Bc", "Bbaryon", "induced_tau", "Kmeson", "WZ"]
else:
    sys.exit(f"unknown grid {grid}")

# combine with the repo code (in the production venv, as published)
code = (f"from grendel.models.hnl.production.combine_channels import combine_for_point as c; "
        f"print(c({flavor!r}, {mass!r}, channels={channels!r}, strict=True))")
r = subprocess.run([PROD, "-c", code], env=env, cwd=job, capture_output=True, text=True)
if r.returncode:
    sys.exit(f"combine failed: {r.stderr[-2000:]}")
comb = vec / flavor / "combined" / f"mN_{lab}.csv"
log(f"combined {r.stdout.strip()} rows")


def digest(p: Path):
    h, n, w = hashlib.sha256(), 0, 0.0
    with open(p, "rb") as f:
        for line in f:
            h.update(line)
            n += 1
            w += float(line.split(b",", 1)[0])
    return {"rows": n, "w_sum": w, "sha256": h.hexdigest()}


summary = {"grid": grid, "flavor": flavor, "mass_GeV": mass, "channels": {}}
for d in sorted((vec / flavor).iterdir()):
    p = d / f"mN_{lab}.csv"
    if d.name != "combined" and p.exists() and p.stat().st_size > 0:
        summary["channels"][d.name] = digest(p)
summary["combined"] = digest(comb)
log("digests")

templates = W / "templates"
tpl = templates / flavor / f"templates_{lab}.npz"
if not tpl.exists():
    sys.exit(f"missing template {tpl}")
summary["template_sha256"] = hashlib.sha256(tpl.read_bytes()).hexdigest()
out_tmp = job / "out"
with open(job / "scan.log", "w") as f:
    r = subprocess.run([sys.executable, "-m", "grendel.models.hnl.scan", "--flavor", flavor, "--mass", mtxt,
                        "--workers", "1", "--thresholds", "3", "10", "--vectors-dir", str(vec),
                        "--templates-dir", str(templates), "--geometry-dir", str(job / "geo"), "--out", str(out_tmp)],
                       env=env, cwd=job, stdout=f, stderr=subprocess.STDOUT)
if r.returncode or not (out_tmp / "sensitivity.csv").exists():
    sys.exit(f"scan failed ({job / 'scan.log'})")
log("scanned: " + (out_tmp / "sensitivity.csv").read_text().splitlines()[-1][:120])
summary["elapsed_s"] = round(time.time() - t0, 1)
tmp = out.with_name(out.name + ".part")
shutil.rmtree(tmp, ignore_errors=True)
tmp.mkdir(parents=True)
for f in ("sensitivity.csv", "run_metadata.json", "scan_status.json"):
    if (out_tmp / f).exists():
        shutil.copy(out_tmp / f, tmp / f)
for f in ("scan.log", "prod_mesons.log"):
    if (job / f).exists():
        shutil.copy(job / f, tmp / f)
(tmp / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
shutil.rmtree(out, ignore_errors=True)
tmp.rename(out)
shutil.rmtree(job, ignore_errors=True)
log(f"-> {out}")
