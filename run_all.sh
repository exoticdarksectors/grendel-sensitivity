#!/usr/bin/env bash
# Everything, after env/bootstrap.sh (and, for SensCalc, the Wolfram Engine install + activation; README.md):
# the five pipelines start together and share the CPU pool of pipelines/common/slot_run.py, then the figures are
# gathered into figures/. Resume-safe: rerun it after an interruption and every pipeline skips its finished work.
#
#   source env/activate.sh
#   setsid nohup bash run_all.sh > logs/run_all.out 2>&1 < /dev/null &
#
# Wall time on a 12-core laptop (Ryzen 5 PRO 7530U, 13 GB): about two days, dominated by the HNL fine grid
# (pipelines/hnl/README.md). Progress: logs/<area>/run_all.{status,log,out}; CPU: pipelines/common/slots_status.sh.
set -uo pipefail
REPRO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$REPRO/env/activate.sh"
cd "$REPRO"
mkdir -p logs/{fonll,scalar,hnl,alp,senscalc}
echo "$(date '+%F %T') run_all start (pid $$)"

start() {  # start <area> <log>
  echo "$(date '+%F %T') start $1 -> $2"
  setsid bash "pipelines/$1/run_all.sh" > "$2" 2>&1 < /dev/null &
  pids+=("$!"); names+=("$1")
}
pids=(); names=()
start fonll  logs/fonll/run_all.log
start scalar logs/scalar/run_all.out
start hnl    logs/hnl/run_all.out
start alp    logs/alp/run_all.out
if [ -x env/wolfram/bin/wolframscript ] && wolframscript -code 'Print[$LicenseID]' > /dev/null 2>&1; then
  start senscalc logs/senscalc/run_all.log
else
  echo "$(date '+%F %T') SKIP senscalc: Wolfram Engine not installed or not activated (pipelines/senscalc/README.md)"
fi

rc=0
for i in "${!pids[@]}"; do
  if wait "${pids[$i]}"; then echo "$(date '+%F %T') done ${names[$i]}"
  else echo "$(date '+%F %T') FAIL ${names[$i]} (see its log)"; rc=1; fi
done
bash figures/collect.sh || rc=1
echo "$(date '+%F %T') run_all end (exit $rc)"
exit $rc
