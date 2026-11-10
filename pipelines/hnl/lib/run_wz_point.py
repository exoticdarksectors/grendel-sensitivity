"""W/Z -> l N at one (flavour, mass): the repo's run_single_point on a copy of the compiled process directory.

    python run_wz_point.py <flavor> <mass> <n_events> <nb_core> <randinit> <vectors_dir> <scratch_dir> [--gzip]

As hnl_fine_condor_20261014/condor/run_wz.py (the fine-grid recipe): get_or_create_process_dir is replaced by a fresh
copy of work/hnl/mg5/master/hnl_<flavor> (step 03), and MadGraph gets up to 6 h. SubProcesses/randinit is set to the
published seed offset before the run (MadGraph's iseed = 0 reads it): r=21 for the June v1 W/Z (coarse grid; the
template value, the directory was rebuilt at every mass), r=27 (Ue, Utau) / r=33 (Umu; Umu 3.845: r=34) for the fine
grid (the copied lxplus directory). nb_core does not change the events (Umu 3 GeV, 1M: --nb_core=1 and 4 give
identical files here), so the pipeline runs every point on one core (the fine grid used 4 on lxplus). Writes <vectors_dir>/<flavor>/WZ/mN_<label>.csv (weights x K = 1.3) and keeps
MadGraph's generate_events log next to it as mN_<label>.generate_events.log. --gzip stores the CSV as .csv.gz
(the fine grid: 1M rows each).
"""
import functools
import gzip
import os
import shutil
import sys
from pathlib import Path

flavor, mass, n_events, nb_core, randinit = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
vectors, scratch = Path(sys.argv[6]).resolve(), Path(sys.argv[7]).resolve()
os.environ["GRENDEL_HNL_VECTORS_DIR"] = str(vectors)

import grendel.models.hnl.production.madgraph.run_wz_production as wz  # noqa: E402
from grendel.models.hnl.production.madgraph.runner import run_events  # noqa: E402
from grendel.models.hnl.production.paths import llp_csv_path  # noqa: E402
from grendel.production.constants import FLAVOR_TO_MG5  # noqa: E402

master = Path(os.environ["REPRO"]) / "work/hnl/mg5/master" / f"hnl_{FLAVOR_TO_MG5[flavor]}"
assert (master / "bin/generate_events").exists() and list(master.glob("SubProcesses/P*/madevent")), master
proc = scratch / master.name
shutil.rmtree(proc, ignore_errors=True)
scratch.mkdir(parents=True, exist_ok=True)
shutil.copytree(master, proc, symlinks=True)
for p in (proc / "Events").glob("*"):
    shutil.rmtree(p) if p.is_dir() else p.unlink()
(proc / "SubProcesses/randinit").write_text(f"r={randinit}\n")

wz.get_or_create_process_dir = lambda fl: proc
wz.mg5_run_events = functools.partial(run_events, timeout=6 * 3600)
ok = wz.run_single_point(flavor, mass, n_events, nb_core=nb_core)
csv = llp_csv_path(flavor, "WZ", mass)
log = proc / f"generate_events_run_{csv.stem.removeprefix('mN_')}.log"
if log.exists():
    shutil.copy(log, csv.with_suffix(".generate_events.log"))
shutil.rmtree(proc, ignore_errors=True)
if ok and "--gzip" in sys.argv[8:]:
    gz = csv.with_name(csv.name + ".gz")
    with open(csv, "rb") as fi, gzip.open(str(gz) + ".part", "wb", compresslevel=1) as fo:
        shutil.copyfileobj(fi, fo)
    os.replace(str(gz) + ".part", gz)
    csv.unlink()
sys.exit(0 if ok else 1)
