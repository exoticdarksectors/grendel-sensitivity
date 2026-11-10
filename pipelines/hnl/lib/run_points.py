"""Run the per-point jobs of one stage through the shared CPU pool, resume-safe.

    python run_points.py wz|scan [--grid coarse fine] [--flavor ...] [--jobs N] [--slots K] [--only F:M ...]
                         [--reverse] [--give-up-when FILE ...]

wz    MadGraph W/Z per (flavour, mass): coarse 100k events, r=21 (June); fine 1M events, r=27/33 (Umu 3.845: 34),
      gzipped. One slot per job (--nb_core=1; the fine grid ran --nb_core=4 on lxplus: same events).
scan  lib/scan_point.py per (flavour, mass), one slot each. A point waits (without a slot) until its inputs exist
      (template; coarse: the channel files of steps 05 and 06; fine: the W/Z sample), so scans can follow
      the W/Z production point by point; a fine point also waits (without a slot) until lib/memgate.py would let
      it start. --give-up-when FILE...: once all these files exist (run_all touches them when the producing
      stages have ended) a point whose inputs are still missing fails instead of waiting for ever.
At most N jobs of this stage wait for or hold slots at once; --reverse takes the points in reverse order (a second
driver next to a running one: the per-point locks keep them off each other's points, they meet in the middle).
Logs: logs/hnl/<stage>/<grid>_<flavor>_<mass>.log.
"""
import argparse
import fcntl
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPRO = Path(os.environ["REPRO"])
sys.path.insert(0, str(REPRO / "pipelines/hnl/lib"))
import grids  # noqa: E402
import memgate  # noqa: E402

W = REPRO / "work/hnl"
LIB = REPRO / "pipelines/hnl/lib"
SLOT = [sys.executable, str(REPRO / "pipelines/common/slot_run.py")]

ap = argparse.ArgumentParser()
ap.add_argument("stage", choices=["wz", "scan"])
ap.add_argument("--grid", nargs="+", default=["coarse", "fine"], choices=["coarse", "fine"])
ap.add_argument("--flavor", nargs="+", default=grids.FLAVORS)
ap.add_argument("--jobs", type=int, default=2)
ap.add_argument("--slots", type=int, default=1)
ap.add_argument("--nb-core", type=int, default=1, help="MadGraph nb_core (needs as many slots)")
ap.add_argument("--only", nargs="+", default=None, help="restrict to flavor:mass points")
ap.add_argument("--reverse", action="store_true", help="points in reverse order")
ap.add_argument("--give-up-when", nargs="+", default=[], type=Path,
                help="scan: fail a point whose inputs are missing once all these files exist")
args = ap.parse_args()

points = [(g, fl, m) for g in args.grid for fl in args.flavor for m in (grids.COARSE if g == "coarse" else grids.FINE)]
if args.only:
    want = {(o.split(":")[0], round(float(o.split(":")[1]), 3)) for o in args.only}
    points = [p for p in points if (p[1], round(p[2], 3)) in want]
if args.reverse:
    points.reverse()


def wz_out(g, fl, m):
    base = W / g / "wz/llp_4vectors" / fl / "WZ" / f"mN_{grids.label(m)}.csv"
    return base.with_name(base.name + ".gz") if g == "fine" else base


def wz_done(g, fl, m):
    """The W/Z output is complete: run_wz_point copies MadGraph's log after the K-factor rescaling (coarse) and
    writes the .gz by rename after it (fine)."""
    log = W / g / "wz/llp_4vectors" / fl / "WZ" / f"mN_{grids.label(m)}.generate_events.log"
    return wz_out(g, fl, m).exists() and log.exists()


def ready(g, fl, m):
    lab = grids.label(m)
    if not wz_done(g, fl, m):
        return False
    tpl = W / "templates" / fl / f"templates_{lab}.npz"
    if not tpl.exists() or time.time() - tpl.stat().st_mtime < 120:      # a bundle being written is not ready
        return False
    need = []
    if g == "coarse":
        need += [W / "coarse" / prod / "llp_4vectors" / fl / ch / f"mN_{lab}.csv"
                 for prod, chans in (("asbuilt", ["Dmeson", "Bc", "Kmeson"]),
                                     ("v1", ["Bmeson", "Bbaryon", "tau", "induced_tau"])) for ch in chans]
        need += [W / "coarse" / prod / f".done_{fl}" for prod in ("asbuilt", "v1")]
    return all(p.exists() for p in need)


def job(p):
    """One point under a per-point lock (work/hnl/locks), so two drivers never run the same point."""
    g, fl, m = p
    locks = W / "locks"
    locks.mkdir(parents=True, exist_ok=True)
    with open(locks / f"{args.stage}_{g}_{fl}_{m:g}.lock", "a+") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        return _job(p)


def _job(p):
    g, fl, m = p
    tag = f"{g}_{fl}_{m:g}"
    if args.stage == "wz":
        if wz_done(g, fl, m):
            return tag, 0, 0.0, "exists"
        nev, r = (100_000, grids.RANDINIT_COARSE) if g == "coarse" else (1_000_000, grids.fine_randinit(fl, m))
        nbc = args.nb_core
        cmd = [sys.executable, str(LIB / "run_wz_point.py"), fl, f"{m:g}", str(nev), str(nbc), str(r),
               str(W / g / "wz/llp_4vectors"), str(W / "mg5/runs" / tag)] + (["--gzip"] if g == "fine" else [])
        k = max(args.slots, nbc)
    else:
        if (W / "results" / g / f"{fl}_{m:g}" / "sensitivity.csv").exists():
            return tag, 0, 0.0, "exists"
        while not ready(g, fl, m):
            if args.give_up_when and all(f.exists() for f in args.give_up_when):
                time.sleep(150)                       # a template written just before the marker needs 120 s
                if not ready(g, fl, m):
                    print(f"{time.strftime('%F %T')} {tag}: inputs missing and their stages have ended", flush=True)
                    return tag, 1, 0.0, "ran"
                break
            time.sleep(60)
        while g == "fine" and not memgate.would_pass():
            time.sleep(60)
        cmd = [sys.executable, str(LIB / "scan_point.py"), g, fl, f"{m:g}"]
        k = 1
    logd = REPRO / "logs/hnl" / args.stage
    logd.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with open(logd / f"{tag}.log", "w") as f:
        rc = subprocess.run(SLOT + ["-n", str(k), "--label", f"hnl:{args.stage}:{tag}", "--"] + cmd,
                            stdout=f, stderr=subprocess.STDOUT).returncode
    return tag, rc, time.time() - t0, "ran"


todo = len(points)
print(f"{args.stage}: {todo} points, {args.jobs} concurrent jobs x {args.slots if args.stage == 'wz' else 1} slot(s)",
      flush=True)
fails = 0
with ThreadPoolExecutor(max_workers=args.jobs) as ex:
    futs = [ex.submit(job, p) for p in points]
    for i, f in enumerate(as_completed(futs), 1):
        tag, rc, dt, what = f.result()
        if what == "ran":
            print(f"{time.strftime('%F %T')} [{i}/{todo}] {tag}: {'OK' if rc == 0 else f'FAIL rc={rc}'} "
                  f"({dt:.0f} s incl. queue)", flush=True)
        fails += rc != 0
print(f"{args.stage}: done, {fails} failures", flush=True)
sys.exit(1 if fails else 0)
