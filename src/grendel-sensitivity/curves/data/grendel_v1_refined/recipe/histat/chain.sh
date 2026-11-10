#!/bin/bash
# 4x production statistics for the HNL tip window 2.6-3.7 GeV (Ue, Umu, Utau):
# B, B_c, Lambda_b pools of 400k and 400k MadGraph W/Z events per mass, seed 4242.
# One scan worker (~7 GB) at a time; a memory guard stops a stage below 12 % free.
set -uo pipefail
H=${HISTAT_DIR:?campaign directory}
REPO=${GRENDEL_REPO:?this repository}
CODE_V1=${CODE_V1:?worktree of b6670a8}
HNLPY=${HNL_CONDA_PYTHON:?python of the hnl conda env}
PY=$HOME/.venvs/grendel-sensitivity/bin/python
MW="2.6 2.8 3 3.05 3.1 3.15 3.2 3.25 3.3 3.35 3.4 3.45 3.5 3.55 3.6 3.62 3.63 3.64 3.65 3.675 3.685 3.7"
STATUS=$H/chain.status
export PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=$H/.mpl GRENDEL_THIRD_PARTY_DIR=$H/third_party
mark() { echo "$(date '+%F %T')  $*" | tee -a $STATUS; }
freepct() { memory_pressure -Q 2>/dev/null | awk -F': ' '/free percentage/ {gsub("%","",$2); print $2+0}'; }
run() {  # run <name> <log> <command...>, guarded
  local name=$1 log=$2; shift 2
  mark "START  $name"
  "$@" >> "$log" 2>&1 &
  local pid=$!
  while kill -0 $pid 2>/dev/null; do
    local f=$(freepct)
    if [ -n "$f" ] && [ "$f" -lt 12 ]; then
      mark "GUARD  $name: memory free ${f} %, stopping"; pkill -TERM -P $pid; kill -TERM $pid; sleep 5; pkill -KILL -P $pid; return 9
    fi
    sleep 5
  done
  wait $pid; local rc=$?
  mark "$( [ $rc -eq 0 ] && echo 'OK    ' || echo 'FAIL  ' )$name (exit=$rc)"
  return $rc
}
cd $REPO
[ -d $H/code_main ] || git -C $REPO worktree add --detach $H/code_main d717b24 >> $H/logs/setup.log 2>&1
mark "chain started (production: main d717b24; v1 scan: b6670a8; v2 scan: d717b24)"

run production $H/logs/production.log env GRENDEL_MG5_EXE=$H/third_party/MG5_aMC_v3_6_6/bin/mg5_aMC \
    GRENDEL_LHAPDF_CONFIG=${GRENDEL_LHAPDF_CONFIG:?} \
    GRENDEL_HNL_VECTORS_DIR=$H/llp_4vectors GRENDEL_WORK_DIR=$H/work PYTHONPATH=$H/code_main PYTHONUNBUFFERED=1 \
    $HNLPY -m grendel.models.hnl.production.run_all --flavor Ue Umu Utau --masses $MW \
    --channels bottom bc baryon wz --no-prompt-tau --skip-combine \
    --n-pool 400000 --wz-nevents 400000 --wz-nb-core 4 --workers 6 --seed 4242 || exit 1

run combine $H/logs/combine.log env GRENDEL_HNL_VECTORS_DIR=$H/llp_4vectors PYTHONPATH=$H/code_main $PY -c "
from grendel.models.hnl.production.combine_channels import combine_for_point
from grendel.models.hnl.production.paths import llp_csv_path
MW = [float(x) for x in '$MW'.split()]
for fl in ('Ue', 'Umu', 'Utau'):
    for m in MW:
        for ch in ('Dmeson', 'induced_tau', 'Kmeson'):   # closed above 2.6 GeV: empty, as in the published production
            p = llp_csv_path(fl, ch, m)
            if not p.exists(): p.write_text('')
        combine_for_point(fl, m, channels=['Bmeson', 'Dmeson', 'Bc', 'Bbaryon', 'induced_tau', 'Kmeson', 'WZ'], strict=True)
print('combined', 3 * len(MW))
" || exit 1

run assemble $H/logs/assemble.log env PYTHONPATH=$H/code_main $PY $H/assemble.py "$MW" || exit 1

run "v1 scan" $H/logs/scan_v1.log env PYTHONPATH=$CODE_V1 sh -c "cd $CODE_V1 && $PY -m grendel.models.hnl.scan \
    --flavor Ue Umu Utau --mass $MW --workers 1 --thresholds 3 10 \
    --vectors-dir $H/vectors_v1 --templates-dir $H/templates_v1 --geometry-dir $H/geo_v1 --out $H/results_v1" || exit 1
run "v2 scan" $H/logs/scan_v2.log env PYTHONPATH=$H/code_main sh -c "cd $H/code_main && $PY -m grendel.models.hnl.scan \
    --flavor Ue Umu Utau --mass $MW --workers 1 --thresholds 3 10 \
    --vectors-dir $H/vectors_v2 --templates-dir $H/templates_v2 --geometry-dir $H/geo_v2 --out $H/results_v2" || exit 1
mark "chain finished"
