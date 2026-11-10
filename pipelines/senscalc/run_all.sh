#!/usr/bin/env bash
# run_all.sh -- the whole Wolfram/SensCalc chain, resume-safe. Prerequisite: 00_install_wolfram.sh, 01_fetch_sources.sh,
# and a one-time `wolframscript -activate` with the user's Wolfram ID.
#
#   cd $REPRO && source env/activate.sh
#   setsid nohup bash pipelines/senscalc/run_all.sh > logs/senscalc/run_all.log 2>&1 < /dev/null &
#
# Every SensCalc run takes $SENSCALC_SUBKERNELS (default 4) slots of the shared pool through slot_run.py -n 4, one run
# at a time (they share one SensCalc tree, as in the campaign); the light steps take 1 slot. Steps already done are skipped (work/senscalc/markers);
# to redo one, delete its marker. SKIP_SMOKE=1 skips 10_smoke.sh.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
: "${REPRO:?source env/activate.sh first}"
echo "run_all start $(date '+%F %T')"
bash 02_check_engine.sh
bash 22_check_quartic.sh
bash 20_export_senscalc_2501.sh
bash 21_export_senscalc_2501_production.sh
bash 23_check_scalar_1809.sh || echo "23_check_scalar_1809.sh reported differences (see work/senscalc/validation)"
[ "${SKIP_SMOKE:-0}" = 1 ] || bash 10_smoke.sh
bash 11_nb1_acceptances.sh
bash 12_nb2_distributions.sh
bash 13_nb3_sensitivity.sh
bash 14_curves_validate.sh
echo "run_all done $(date '+%F %T')"
