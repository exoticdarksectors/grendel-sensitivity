#!/usr/bin/env bash
# The whole BC4/BC5 chain, resume-safe: every step skips the work already done, so rerunning this script after an
# interruption continues where it stopped. Status lines go to logs/scalar/run_all.status.
#
#   source env/activate.sh
#   setsid nohup bash pipelines/scalar/run_all.sh > logs/scalar/run_all.out 2>&1 &
#
# Order: inputs (00 setup, 01 BC4 four-vectors, 05 BC5 four-vectors, 06 exHad templates), then the scans (07 control
# masses, 08 every mass; SCALAR_JOBS jobs in flight, one slot each), while the Higgs-histogram chain (02 MG5, 03
# shower, 04 histogram), which validates the tracked histogram the BC5 four-vectors use, runs alongside; then
# 09 tables + comparisons and 10 figures (they need only 08), and last 11 (BC5 rows with the rebuilt vs the tracked
# Higgs histogram, which needs the Higgs chain). Knobs: SCALAR_JOBS (scan jobs in flight, default 4 in 07 and 6 in
# 08; each holds one slot of the shared pool), SCALAR_MG5_SLOTS (default 4), SCALAR_TPL_SLOTS (default 4).
# A failing step stops the script (the status file says which); rerunning it resumes.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
ST="$L/run_all.status"
step() {
  echo "$(date '+%F %T') START $1" >> "$ST"
  if bash "$S/$1"; then echo "$(date '+%F %T') OK    $1" >> "$ST"; else echo "$(date '+%F %T') FAIL  $1" >> "$ST"; return 1; fi
}
step 00_setup.sh
step 01_bc4_vectors.sh
step 05_bc5_vectors.sh
step 06_templates.sh
( step 02_bc5_higgs_mg5.sh && step 03_bc5_higgs_shower.sh && step 04_bc5_higgs_hist.sh ) > "$L/run_all_higgs.out" 2>&1 &
higgs=$!
step 07_scan_controls.sh
step 08_scan_all.sh
step 09_collect.sh
step 10_figures.sh
wait "$higgs" || true
step 11_higgs_hist_impact.sh
echo "$(date '+%F %T') run_all finished" >> "$ST"
