#!/usr/bin/env bash
# Step 04: the central meson grids from the raw quark grids of step 03, with the shipped generator
#     python -m grendel.production.fonll_grids.generate --pdf nlo --quark <q> --reuse-existing-grids
# which finds the complete raw grid of step 03 (gen.raw_grid_complete) and runs only fragmfonll and the writer:
# bottom = Kartvelishvili alpha = 24.2 (mode 2); charm = BCFY pseudoscalar r = 0.1 (mode 5) + D* vector (mode 4)
# at the D* parent kinematics, combined with the step 02 feeddown weight. Minutes per quark, one slot each.
#
#   source env/activate.sh && bash pipelines/fonll/04_fragment_central.sh [bottom] [charm]
#
# Output: work/fonll/grids/output/fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_{bottom,charm}.dat
# (+ generation_summary_nlo_<q>.json), the files GRENDEL_FONLL_{BOTTOM,CHARM}_GRID would point at.
set -euo pipefail
: "${REPRO:?source env/activate.sh first}"
QUARKS=("$@"); [ ${#QUARKS[@]} -gt 0 ] || QUARKS=(bottom charm)
LOG="$REPRO/logs/fonll/04_fragment_central.log"
mkdir -p "$(dirname "$LOG")"
exec > >(tee -a "$LOG") 2>&1
echo "=== $(date -Is) 04_fragment_central: ${QUARKS[*]}"
T="$REPRO/pipelines/fonll/task.sh"
pids=(); rc=0
for q in "${QUARKS[@]}"; do bash "$T" "frag:$q" & pids+=($!); done
for p in "${pids[@]}"; do wait "$p" || rc=1; done
echo "=== $(date -Is) 04_fragment_central done (rc=$rc)"
exit $rc
