#!/usr/bin/env bash
# One unit of FONLL work, used by steps 02-04 and run_all.sh. Each unit holds one CPU slot (slot_run -n 1) and a
# per-unit flock (work/fonll/locks/<unit>.lock), so two drivers never compute the same unit twice: the second
# waits, then finds the unit complete and returns. Every unit is resume-safe on its own.
#
#   task.sh calib              charm CTEQ6.6 D* feeddown calibration (generate --calibrate-only)
#   task.sh closure            bottom CTEQ6.6 public-point closure (closure.py)
#   task.sh raw:<quark>:<i>    rapidity chunk i of the NNPDF4.0 raw quark grid (raw_chunk.py, FONLL_NCHUNKS chunks)
#   task.sh merge:<quark>      merge the chunks into run/nlo_<quark>_central/nl_<q>_.out (no slot: seconds)
#   task.sh frag:<quark>       generate --pdf nlo --quark <quark> --reuse-existing-grids (fragmentation + grid file)
set -euo pipefail
: "${REPRO:?source env/activate.sh first}"
source "$REPRO/pipelines/fonll/env.sh"
cd "$REPRO"
unit="$1"
NCH="${FONLL_NCHUNKS:-4}"
SLOT=(python "$REPRO/pipelines/common/slot_run.py" -n 1 --label "fonll:$unit" --)
WS=(--fonll "$GRENDEL_FONLL_DIR" --out-dir "$GRENDEL_FONLL_GRID_DIR" --run-dir "$FONLL_RUN_DIR" --log-dir "$FONLL_LOG_DIR")
mkdir -p "$FONLL_WORK/locks" "$GRENDEL_FONLL_GRID_DIR" "$FONLL_RUN_DIR" "$FONLL_LOG_DIR" "$FONLL_WORK/validation"
exec 9>"$FONLL_WORK/locks/${unit//:/_}.lock"
flock 9
t0=$SECONDS
stamp() { echo "[$(date -Is)] $unit: $*"; }

grid_done() {
  local f="$GRENDEL_FONLL_GRID_DIR/fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_$1.dat"
  [ -f "$f" ] && [ "$(grep -vc '^#' "$f")" -eq 10000 ]
}

case "$unit" in
  calib)
    json="$GRENDEL_FONLL_GRID_DIR/charm_feeddown_calibration.json"
    if [ -f "$json" ]; then stamp "done already"; exit 0; fi
    raw="$FONLL_RUN_DIR/validate_cteq66_charm_y0/val_c_.out"
    # ensure_charm_calibration_grid reuses any existing val_c_.out: drop one an interrupted run left short
    if [ -f "$raw" ] && [ "$(wc -l < "$raw")" -ne 241 ]; then stamp "removing incomplete $raw"; rm -f "$raw"; fi
    stamp "start"
    "${SLOT[@]}" python -m grendel.production.fonll_grids.generate "${WS[@]}" --calibrate-only
    ;;
  closure)
    json="$FONLL_WORK/validation/bottom_cteq66_public_point.json"
    if [ -f "$json" ]; then stamp "done already"; exit 0; fi
    stamp "start"
    "${SLOT[@]}" python "$REPRO/pipelines/fonll/closure.py" --out "$json"
    ;;
  raw:*)
    IFS=: read -r _ q i <<< "$unit"
    if grid_done "$q"; then stamp "grid exists"; exit 0; fi
    if python "$REPRO/pipelines/fonll/raw_chunk.py" --quark "$q" --nchunks "$NCH" --chunk "$i" --check-only 2>/dev/null; then
      stamp "done already"; exit 0
    fi
    stamp "start"
    "${SLOT[@]}" python "$REPRO/pipelines/fonll/raw_chunk.py" --quark "$q" --nchunks "$NCH" --chunk "$i"
    ;;
  merge:*)
    q="${unit#merge:}"
    if grid_done "$q"; then stamp "grid exists"; exit 0; fi
    python "$REPRO/pipelines/fonll/raw_chunk.py" --quark "$q" --nchunks "$NCH" --merge
    ;;
  frag:*)
    q="${unit#frag:}"
    if grid_done "$q"; then stamp "grid exists"; exit 0; fi
    python "$REPRO/pipelines/fonll/raw_chunk.py" --quark "$q" --nchunks "$NCH" --merge
    if [ "$q" = charm ] && [ ! -f "$GRENDEL_FONLL_GRID_DIR/charm_feeddown_calibration.json" ]; then
      stamp "no charm_feeddown_calibration.json (task calib)"; exit 1
    fi
    stamp "start"
    "${SLOT[@]}" python -m grendel.production.fonll_grids.generate "${WS[@]}" --pdf nlo --quark "$q" --reuse-existing-grids
    grid_done "$q"
    ;;
  *) echo "unknown unit $unit" >&2; exit 2 ;;
esac
stamp "finished in $((SECONDS - t0)) s"
