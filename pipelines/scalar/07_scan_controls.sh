#!/usr/bin/env bash
# Step 07 -- control masses at full statistics, before the production (08).
#
#   BC4 0.5 1 2 3.6 (1.0 = the deepest point, 3.6 = the reach), BC5 0.3 1 15 50 (50 = the deepest point)
# each scanned twice with the published settings (scan_one.sh):
#   control_reftpl -- with the PUBLISHED exHad bundles (macOS arm64): only the four-vectors are ours, so the row must
#                     equal the published one to the precision of the vectors (libm last-digit noise, ~1e-8);
#                     this checks the four-vectors (incl. the BC4 particle-0.26.2 weights) and the scan on Linux;
#   main           -- with our Linux exHad templates (06): the production rows of these masses; differ from the
#                     published ones by the decay Monte Carlo only.
# Report: work/scalar/validation/controls_<bench>_<run>.txt (+ .tsv). Needs 01, 05 and 06. control_reftpl needs the
# published bundles of the reference data ($GRENDEL_REF) and is skipped without them.
# Up to SCALAR_JOBS (default 4) jobs at once, one slot each; ~15-25 min per job.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
step_lock
J=${SCALAR_JOBS:-4}
CTRL_BC4="0.5 1 2 3.6"; CTRL_BC5="0.3 1 15 50"
RUNS="main"
ref_ok "control_reftpl scans (the published exHad bundles)" "$REF_TEMPLATES" && RUNS="control_reftpl main"
{
  for r in $RUNS; do
    for m in $CTRL_BC4; do echo "bc4 $m $r"; done
    for m in $CTRL_BC5; do echo "bc5 $m $r"; done
  done
} | xargs -P "$J" -L 1 bash "$S/scan_one.sh" 2>&1 | tee -a "$L/07_scan_controls.log"
mkdir -p "$W/validation"
for b in bc4 bc5; do
  for run in $RUNS; do
    python "$S/py/compare_rows.py" "$b" "$PUB/$b/sensitivity.csv" "$W/$b/scans/$run/out" \
      --tsv "$W/validation/controls_${b}_${run}.tsv" | tee "$W/validation/controls_${b}_${run}.txt"
  done
done
step_ok 07_scan_controls
