#!/usr/bin/env bash
# Step 02: the CTEQ6.6 inputs (two small FONLL runs, ~15-25 min each, one slot each, run in parallel).
#
#   source env/activate.sh && bash pipelines/fonll/02_calibrate.sh
#
# a) task calib: generate --calibrate-only. FONLL CTEQ6.6 (LHAPDF 10550) charm quark grid at y = -0.5, 0, 0.5
#    (work/fonll/grids/run/validate_cteq66_charm_y0/val_c_.out), the direct BCFY pseudoscalar (mode 5) and D*
#    vector (mode 4) fragmentation at the nine public FONLL D0 points, and the ternary-search fit of the D* -> D0
#    feeddown weight, cached in work/fonll/grids/output/charm_feeddown_calibration.json (read by step 04 charm).
#    Reference: fonll-local/output/charm_feeddown_calibration.json (weight 1.2152533641182024; the tracked charm
#    grid header records 1.21525336412) and fonll-local/run/validate_cteq66_charm_y0/.
# b) task closure (validation only): CTEQ6.6 bottom quark grid at y = -0.5, 0, 0.5 fragmented at pT = 5 GeV,
#    y = 0 (Kartvelishvili alpha = 24.2) -> work/fonll/validation/bottom_cteq66_public_point.json; against the
#    public FONLL web value 7.5329e6 pb/GeV and the original local value 7.535277e6 (data/fonll/PROVENANCE.md).
#    The original grid (fonll-local/run/cteq66_bottom) had 13 rapidities; with ifrframe = 1 the fragmentation at
#    y = 0 reads only the y = 0 row, so 3 rows suffice.
set -euo pipefail
: "${REPRO:?source env/activate.sh first}"
LOG="$REPRO/logs/fonll/02_calibrate.log"
mkdir -p "$(dirname "$LOG")"
exec > >(tee -a "$LOG") 2>&1
echo "=== $(date -Is) 02_calibrate"
T="$REPRO/pipelines/fonll/task.sh"
rc=0
bash "$T" calib & p1=$!
bash "$T" closure & p2=$!
wait $p1 || rc=1
wait $p2 || rc=1
echo "=== $(date -Is) 02_calibrate done (rc=$rc)"
exit $rc
