#!/usr/bin/env bash
# Unattended full FONLL uncertainty campaign.
#
#   Phase 1  alpha_s companion central grids (as_01170 / as_01190 x bottom/charm)
#            through the plain path: their own generation_summary files, so they
#            cannot clobber the campaign manifest. No grid reuse, so a partial raw
#            grid is rebuilt cleanly.
#   Phase 2  alpha_s band -> <grid dir>/envelopes/alphas_manifest.json
#   Phase 3  NNPDF replica campaign in cumulative stages. Each stage rewrites a
#            COMPLETE variation_manifest.json (central + 7 scale + replicas 1..N +
#            2 mass per quark) and refreshes the envelopes, so a coherent set is
#            banked at every stage boundary. Re-fragmenting earlier members is
#            cheap; the cost is the fresh raw grids of the new members.
#
# Prerequisites: the patched FONLL executables (fonll_grids.install --build) and
# the NNPDF4.0 LHAPDF sets incl. as_01170 / as_01190. GRENDEL_FONLL_DIR and
# GRENDEL_FONLL_GRID_DIR (or --out-dir) locate the tree and the output.
#
# Tunables (env): MAXPAR (concurrent grids, default 12), REUSE (default 1),
#   STAGES (cumulative member caps, default "25 40 55 70 85 100"),
#   PY (python with this package installed, default python3).
set -u
PY=${PY:-python3}
MAXPAR=${MAXPAR:-12}
REUSE=${REUSE:-1}
STAGES=${STAGES:-"25 40 55 70 85 100"}
GRID_DIR=${GRENDEL_FONLL_GRID_DIR:-${GRENDEL_WORK_DIR:-grendel_work}/fonll_grids/output}
LOG_DIR=$(dirname "$GRID_DIR")/logs
mkdir -p "$LOG_DIR"
LOG=$LOG_DIR/campaign_$(date +%Y%m%d_%H%M%S).log
exec >>"$LOG" 2>&1

GEN="$PY -m grendel.production.fonll_grids.generate --out-dir $GRID_DIR"
COMBINE="$PY -m grendel.production.fonll_grids.combine"
reuse_flag=""
[ "$REUSE" = "1" ] && reuse_flag="--reuse-existing-grids"

echo "================ campaign start $(date) ================"
echo "grid_dir=$GRID_DIR max_parallel=$MAXPAR reuse='$reuse_flag' stages='$STAGES' log=$LOG"

echo "---- Phase 1: alpha_s grids $(date) ----"
$GEN --pdf nlo_as_01170 --quark bottom --quark charm --grid-workers 6 &
A1=$!
$GEN --pdf nlo_as_01190 --quark bottom --quark charm --grid-workers 6 &
A2=$!
wait $A1; r1=$?
wait $A2; r2=$?
echo "---- alpha_s grids done $(date) (exit $r1/$r2) ----"

echo "---- Phase 2: alpha_s combiner $(date) ----"
$COMBINE alphas --out-dir "$GRID_DIR"

for TARGET in $STAGES; do
  echo "======== Phase 3 stage: replicas 1..${TARGET} $(date) ========"
  $GEN --campaign --scale-variations --mass-variations --pdf-members 1-${TARGET} \
       $reuse_flag --max-parallel "$MAXPAR" --grid-workers 1 --compress-logs
  echo "---- stage ${TARGET} grids done; rebuilding envelopes $(date) ----"
  $COMBINE variations --out-dir "$GRID_DIR"
  echo "======== stage ${TARGET} BANKED (complete manifest through replica ${TARGET}) $(date) ========"
done

echo "================ campaign FINISHED $(date) ================"
