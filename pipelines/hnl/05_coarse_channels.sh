#!/bin/bash
# HNL step 05: the non-MadGraph channels of the coarse grid (the 104 masses per flavour of the lxplus scan),
# exactly as the published inputs were made. Two productions, each over the mass list it was run on (one random
# stream per (flavour, channel) job runs through the masses in order, so the list matters):
#
#  (a) as-built (asbuilt_20261014/hnl/production.sh): grendel run_all with HNLCalc b1d36c6, the repo MASS_GRID
#      (123 masses): charm + B_c (+ kaon, overwritten) at 100k, seed 42, then kaon at 300k, seed 42.
#      Kept: Dmeson, Bc, Kmeson.                                           -> work/hnl/coarse/asbuilt/llp_4vectors
#  (b) June v1 (signal-models/hnl/llpatcolliders_HNL_PR1/hnl/tmp/runs/central_newgrids_20260623/production.log):
#      the code that ran, llpatcolliders 06e5a18 hnl/ (src/llpatcolliders-hnl-06e5a18: HNLCalc e292cef9 vendored,
#      its FONLL grids, SIGMA_BC_PB = 0.9e6 in the induced-tau B_c -> tau nu source), all channels but W/Z,
#      116 masses, 100k, seed 42, with the prompt-tau pool of step 04.
#      Kept: Bmeson, Bbaryon, tau, induced_tau.                            -> work/hnl/coarse/v1/llp_4vectors
#
# run_all gives job i (flavour-major, channels in a fixed order) the seed <seed>+i, so one invocation per flavour
# with the seed its jobs had in the three-flavour command is the same production, split over single-slot jobs:
#   (a) charm/bc/kaon: Ue 42, Umu 45, Utau 48; kaon 300k: Ue 42, Umu 43, Utau 44
#   (b) bottom/charm/bc/induced_tau/kaon/baryon/prompt_tau: Ue 42, Umu 49, Utau 56
# W/Z is step 06. Channel files outside the 104 coarse masses are deleted afterwards (2.0-4.2 GeV is the fine
# grid). Resume-safe: each (production, flavour) writes a .done stamp.   05_coarse_channels.sh [all|asbuilt|v1]
set -euo pipefail
WHAT=${1:-all}
source "$(dirname "$0")/../../env/activate.sh"
C=$REPRO/work/hnl/coarse; mkdir -p $C/asbuilt $C/v1/tmp
PY=$REPRO/work/hnl/venv_prod/bin/python
SLOT="python $REPRO/pipelines/common/slot_run.py -n 1"
L=$REPRO/logs/hnl/05; mkdir -p $L
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 MPLCONFIGDIR=$C/.mpl
J=$REPRO/src/llpatcolliders-hnl-06e5a18/hnl
M116=$($PY -c "import sys; sys.path.insert(0, '$REPRO/pipelines/hnl/lib'); import grids; print(*grids.JUNE116)")

asbuilt() {  # asbuilt <flavor> <seed charm/bc/kaon> <seed kaon 300k>   (flock: one producer per flavour)
  exec {fd}>>$C/asbuilt/.lock_$1; flock $fd
  [ -e $C/asbuilt/.done_$1 ] && return 0
  ( cd $C/asbuilt
    export GRENDEL_HNL_VECTORS_DIR=$C/asbuilt/llp_4vectors GRENDEL_HNL_CACHE_DIR=$C/asbuilt/cache
    $SLOT --label hnl:05a_asbuilt_$1 -- $PY -m grendel.models.hnl.production.run_all --flavor $1 \
      --channels charm kaon bc --no-wz --no-prompt-tau --skip-combine --n-pool 100000 --workers 1 --seed $2 \
      > $L/asbuilt_$1.log 2>&1
    $SLOT --label hnl:05a_asbuilt_kaon_$1 -- $PY -m grendel.models.hnl.production.run_all --flavor $1 \
      --channels kaon --no-wz --no-prompt-tau --skip-combine --n-pool 300000 --workers 1 --seed $3 \
      > $L/asbuilt_kaon_$1.log 2>&1 ) && touch $C/asbuilt/.done_$1
}
v1() {  # v1 <flavor> <seed>
  exec {fd}>>$C/v1/.lock_$1; flock $fd
  [ -e $C/v1/.done_$1 ] && return 0
  [ -s $C/tau_pool.csv ] || { echo "missing $C/tau_pool.csv (step 04)"; return 1; }
  ( cd $J && HNL_LLP_VECTORS_DIR=$C/v1/llp_4vectors HNL_TMP_DIR=$C/v1/tmp HNL_CACHE_DIR=$C/v1/cache \
      HNL_TAU_POOL_CSV=$C/tau_pool.csv $SLOT --label hnl:05b_v1_$1 -- $PY run_all.py --flavor $1 \
      --masses $M116 --n-pool 100000 --workers 1 --seed $2 --no-wz --skip-combine > $L/v1_$1.log 2>&1 ) \
    && touch $C/v1/.done_$1
}
pids=()
if [ $WHAT != v1 ]; then asbuilt Ue 42 42 & pids+=($!); asbuilt Umu 45 43 & pids+=($!); asbuilt Utau 48 44 & pids+=($!); fi
if [ $WHAT != asbuilt ]; then v1 Ue 42 & pids+=($!); v1 Umu 49 & pids+=($!); v1 Utau 56 & pids+=($!); fi
rc=0; for p in "${pids[@]}"; do wait $p || rc=1; done
[ $rc -eq 0 ] || { echo "a production failed: see $L"; exit 1; }
for s in asbuilt/.done_Ue asbuilt/.done_Umu asbuilt/.done_Utau v1/.done_Ue v1/.done_Umu v1/.done_Utau; do
  [ -e $C/$s ] || { echo "$WHAT done; waiting for the other production ($s) before the clean-up"; exit 0; }
done

# keep only the coarse masses and the channels each production contributes
$PY - <<'PY'
import os, sys, shutil
from pathlib import Path
sys.path.insert(0, os.environ["REPRO"] + "/pipelines/hnl/lib")
import grids
C = Path(os.environ["REPRO"]) / "work/hnl/coarse"
keep = {"asbuilt": ["Dmeson", "Bc", "Kmeson"], "v1": ["Bmeson", "Bbaryon", "tau", "induced_tau"]}
labels = {f"mN_{grids.label(m)}.csv" for m in grids.COARSE}
n = 0
for prod, chans in keep.items():
    for fl in grids.FLAVORS:
        base = C / prod / "llp_4vectors" / fl
        for d in list(base.iterdir()):
            if d.name not in chans:
                shutil.rmtree(d); continue
            for f in d.glob("mN_*.csv"):
                if f.name not in labels:
                    f.unlink(); n += 1
        missing = [(ch, l) for ch in chans for l in labels if not (base / ch / l).exists()]
        if missing:
            sys.exit(f"{prod} {fl}: {len(missing)} channel files missing, e.g. {missing[:3]}")
print(f"coarse channel files complete; {n} files outside the coarse grid removed")
PY
