#!/usr/bin/env bash
# The whole FONLL area, in order; resume-safe (every unit skips what is finished; per-unit flocks).
#
#   source env/activate.sh
#   setsid nohup bash pipelines/fonll/run_all.sh > logs/fonll/run_all.log 2>&1 < /dev/null &
#
# 00, 01, then ONE queue of 1-slot units with at most FONLL_PARALLEL (default 4) running: the two step-02 units
# (calib, closure) and the 2 x FONLL_NCHUNKS (default 8) raw-grid chunks of step 03; then the merges, step 04
# (bottom and charm in parallel), the fragmentation cross-check on the stored macOS raw grids (xcheck_fragment.py,
# once; needs the optional reference tree) and step 05. Peak: 4 slots.
set -euo pipefail
: "${REPRO:?source env/activate.sh first}"
cd "$REPRO"
P="$REPRO/pipelines/fonll"
NCH="${FONLL_NCHUNKS:-4}"
echo "=== $(date -Is) run_all (pid $$)"
bash "$P/00_fetch_pdfs.sh"
bash "$P/01_build_fonll.sh"
{ echo calib; echo closure
  for q in bottom charm; do for ((i = 0; i < NCH; i++)); do echo "raw:$q:$i"; done; done; } \
  | xargs -P "${FONLL_PARALLEL:-4}" -n 1 bash "$P/task.sh"
bash "$P/task.sh" merge:bottom
bash "$P/task.sh" merge:charm
bash "$P/04_fragment_central.sh"
source "$P/env.sh"
if [ ! -f "$FONLL_WORK/xcheck/output/fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_charm.dat" ] \
   && have_ref "FONLL xcheck: the stored macOS raw grids fragmented here" "$FONLL_ORIG_WS/run" "$FONLL_ORIG_WS/output"; then
  python "$REPRO/pipelines/common/slot_run.py" -n 1 --label fonll:xcheck -- python "$P/xcheck_fragment.py" \
    > "$REPRO/logs/fonll/xcheck_fragment.log" 2>&1
fi
python "$P/05_compare.py"
echo "=== $(date -Is) run_all done"
