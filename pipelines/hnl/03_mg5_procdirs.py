"""HNL step 03: the MadGraph process directories, built once and compiled, then copied for every run.

    W/Z -> l N   work/hnl/mg5/master/hnl_{electron,muon,tau}   (HeavyN UFO, proc_card_<flavor>.dat)
    prompt tau   work/hnl/mg5/master/tau_pool                   (sm, proc_card_tau_production.dat)

Built with the repo helper grendel.models.hnl.production.madgraph.runner.ensure_process_dir (MG5 3.6.6, the
subprocesses compiled with `make madevent MATRIX=matrix*_orig.o`, as in every published run). The published runs
reused one compiled directory per flavour in the same way: the fine grid copied a precompiled directory into every
HTCondor job (hnl_fine_condor_20261014/condor/job.sh); the June v1 runs rebuilt it at every mass (a repo bug, fixed
in c621279), which is the same directory with the MadGraph template seed r=21. The run steps write
SubProcesses/randinit in their copy, so the seed is the published one whatever the directory has seen.

The first import of the HeavyN model converts it for Python 3 (MadGraph writes the converted model and its
py3_model.pkl cache); that happens here, once, before any parallel run.  Resume-safe: a compiled directory is kept.
"""
import os
import sys
from pathlib import Path

REPRO = Path(os.environ["REPRO"])
MG5 = REPRO / "work/hnl/mg5"
os.environ["GRENDEL_HNL_MG5_WORK_DIR"] = str(MG5 / "master")

from grendel.models.hnl.production.madgraph import run_wz_production as wz  # noqa: E402
from grendel.models.hnl.production.madgraph import run_tau_production as tp  # noqa: E402

ok = True
for label, build in [("electron", lambda: wz.get_or_create_process_dir("Ue")),
                     ("muon", lambda: wz.get_or_create_process_dir("Umu")),
                     ("tau", lambda: wz.get_or_create_process_dir("Utau")),
                     ("tau_pool", tp._build_process_dir)]:
    d = build()
    built = d is not None and list(d.glob("SubProcesses/P*/madevent"))
    print(f"{label:9s} {d} subprocesses with madevent: {len(built) if built else 0}", flush=True)
    ok &= bool(built)
sys.exit(0 if ok else 1)
