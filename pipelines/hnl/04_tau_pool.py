"""HNL step 04: the prompt-tau pool (stage 1 of the prompt-tau channel) of the June v1 production.

MadGraph 3.6.6, sm, p p > W -> tau nu and p p > tau tau (proc_card_tau_production.dat), 100 000 events
(N_EVENTS_DEFAULT), seed offset r=21 (the June pool's generate_events_run_pool.log: "Using random number seed
offset = 21", a fresh process directory), K = 1.3 per row: the repo's generate_tau_pool on a copy of
work/hnl/mg5/master/tau_pool. Writes work/hnl/coarse/tau_pool.csv (+ the MadGraph log). Resume-safe.

    python 04_tau_pool.py [nb_core]
"""
import os
import shutil
import sys
from pathlib import Path

REPRO = Path(os.environ["REPRO"])
OUT = REPRO / "work/hnl/coarse"
POOL = OUT / "tau_pool.csv"
SCRATCH = REPRO / "work/hnl/mg5/run_tau_pool"
nb_core = int(sys.argv[1]) if len(sys.argv) > 1 else 1
if POOL.exists() and POOL.stat().st_size > 0:
    print(f"exists: {POOL}")
    sys.exit(0)
OUT.mkdir(parents=True, exist_ok=True)
shutil.rmtree(SCRATCH, ignore_errors=True)
SCRATCH.mkdir(parents=True)
proc = SCRATCH / "tau_pool"
shutil.copytree(REPRO / "work/hnl/mg5/master/tau_pool", proc, symlinks=True)
for p in (proc / "Events").glob("*"):
    shutil.rmtree(p) if p.is_dir() else p.unlink()
(proc / "SubProcesses/randinit").write_text("r=21\n")
os.environ["GRENDEL_HNL_MG5_WORK_DIR"] = str(SCRATCH)
tmp = OUT / "tau_pool.csv.part"
os.environ["GRENDEL_HNL_TAU_POOL_CSV"] = str(tmp)

from grendel.models.hnl.mass_grid import N_EVENTS_DEFAULT  # noqa: E402
from grendel.models.hnl.production.madgraph import run_tau_production as tp  # noqa: E402

ok = tp.generate_tau_pool(N_EVENTS_DEFAULT, nb_core=nb_core, pool_path=tmp)
log = proc / "generate_events_run_pool.log"
if log.exists():
    shutil.copy(log, OUT / "tau_pool.generate_events.log")
if not ok:
    sys.exit(1)
tmp.rename(POOL)
shutil.rmtree(SCRATCH, ignore_errors=True)
print(f"wrote {POOL}")
