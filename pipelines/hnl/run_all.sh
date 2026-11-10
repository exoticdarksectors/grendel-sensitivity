#!/bin/bash
# HNL (BC6-BC8, Fig. 6): the whole chain, resume-safe and safe to start twice (per-point and per-production locks).
# Every CPU-heavy command queues through pipelines/common/slot_run.py. Needs env/bootstrap.sh (toolchain, MadGraph,
# Pythia, exHad, HNLimits). Step 10 (validation) also compares with the campaign outputs of the optional reference
# tree $GRENDEL_REF where available, and reports those comparisons as skipped otherwise (reference/README.md).
#
#   setsid nohup bash pipelines/hnl/run_all.sh > logs/hnl/run_all.out 2>&1 &      status: logs/hnl/run_all.status
#
# Branches, in parallel (concurrent jobs per branch: HNL_JOBS_WZ_COARSE, HNL_JOBS_WZ_FINE, HNL_JOBS_SCAN_COARSE,
# HNL_JOBS_SCAN_FINE; a waiting scan holds no slot; fine scans also pass lib/memgate.py, 3.2 GB each):
#   A  04 tau pool -> 05 coarse channels -> 06 coarse W/Z (312 x 100k)
#   B  07 templates (6 jobs)
#   C  06 fine W/Z (411 x 1M)
#   D  08 coarse scans, E 08 fine scans (each point starts when its inputs exist; once the branches producing them
#      have ended, a point whose inputs are still missing fails instead of waiting: markers in logs/hnl/.ended)
# then 09 table, 10 validation, 11 Fig. 6. Wall time on the 12-core laptop: dominated by the fine grid (411 x ~25 min
# MadGraph + 411 x ~13 min scans, one slot each), i.e. ~260 slot-hours of the ~310 in total.
set -uo pipefail
H=$(cd "$(dirname "$0")" && pwd)
source $H/../../env/activate.sh
L=$REPRO/logs/hnl
S=$L/run_all.status
mark() { echo "$(date '+%F %T')  $*" | tee -a $S; }
step() { local name=$1; shift; mark "START $name"; "$@"; local rc=$?; mark "$([ $rc -eq 0 ] && echo OK || echo FAIL) $name (exit $rc)"; return $rc; }
SLOT="python $REPRO/pipelines/common/slot_run.py -n 1"

mark "run_all started (pid $$)"
step 00_env bash $H/00_env.sh >> $L/00_env.log 2>&1 || exit 1
step 01_kaon_spectrum bash $H/01_kaon_spectrum.sh >> $L/01_kaon_spectrum.out 2>&1 || exit 1
step 02_export_v1_code bash $H/02_export_v1_code.sh >> $L/02_export_v1_code.log 2>&1 || exit 1
step 03_mg5_procdirs $SLOT --label hnl:03_mg5_procdirs -- python $H/03_mg5_procdirs.py >> $L/03_mg5_procdirs.log 2>&1 || exit 1

E=$L/.ended; rm -rf $E; mkdir -p $E
( step 04_tau_pool $SLOT --label hnl:04_tau_pool -- python $H/04_tau_pool.py 1 >> $L/04_tau_pool.log 2>&1 \
  && step 05_coarse_channels bash $H/05_coarse_channels.sh >> $L/05_coarse_channels.out 2>&1 \
  && step 06_wz_coarse bash $H/06_wz.sh coarse ${HNL_JOBS_WZ_COARSE:-2} >> $L/06_wz_coarse.out 2>&1
  rc=$?; touch $E/coarse_inputs; exit $rc ) & A=$!
( step 07_templates bash $H/07_templates.sh >> $L/07_templates.out 2>&1; rc=$?; touch $E/templates; exit $rc ) & B=$!
( step 06_wz_fine bash $H/06_wz.sh fine ${HNL_JOBS_WZ_FINE:-4} >> $L/06_wz_fine.out 2>&1; rc=$?; touch $E/fine_inputs
  exit $rc ) & C=$!
( step 08_scans_coarse bash $H/08_scans.sh coarse ${HNL_JOBS_SCAN_COARSE:-2} --give-up-when $E/coarse_inputs $E/templates \
  >> $L/08_scans_coarse.out 2>&1 ) & D=$!
( step 08_scans_fine bash $H/08_scans.sh fine ${HNL_JOBS_SCAN_FINE:-2} --give-up-when $E/fine_inputs $E/templates \
  >> $L/08_scans_fine.out 2>&1 ) & F=$!
for p in $A $B $C $D $F; do wait $p; done
step 09_collect python $H/09_collect.py >> $L/09_collect.out 2>&1 || { mark "table incomplete: rerun run_all.sh"; exit 1; }
step 10_validate $SLOT --label hnl:10_validate -- python $H/10_validate.py > $L/10_validate.out 2>&1
step 11_plot bash $H/11_plot.sh both >> $L/11_plot.out 2>&1
mark "run_all done"
