#!/usr/bin/env bash
# Step 04 -- the BC10 sensitivity scan, one job per mass, one worker each (bc10_v24_20261019/lxplus/job_bc10.sh):
#   python -m grendel.models.alp.scan --mass M --workers 1 --thresholds 3 10 \
#          --vectors-dir <llp_4vectors> --templates-dir <templates> --geometry-dir <geo> --out <out>
#
#   04_scan.sh VARIANT [J] [mass ...]
#     VARIANT  main     the templates of step 03 (work/alp/templates_exhad)      -> work/alp/scan_main/bc10_<M>/
#              control  the stored macOS templates of the published run ($GRENDEL_REF, optional: without them the
#                       variant prints SKIP and exits 0) -> work/alp/scan_control/bc10_<M>/   (isolates the
#                       four-vector and platform effects)
#     J        jobs in flight (each takes one slot through slot_run; default 4)
#     masses   default: the 153 masses of points_all.txt
#
# Resume: a mass whose bc10_<M>/sensitivity.csv exists is skipped; a job writes into bc10_<M>.tmp and renames it
# when the scan succeeded, so an interrupted job leaves no row. ~6-8 min and <1 GB RAM per mass (one core).
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
VARIANT="${1:?main|control}"; J="${2:-4}"; shift 2 || shift $#
case "$VARIANT" in
  main)    TPL="$W/templates_exhad" ;;
  control) TPL="$CAMP/templates_exhad"
           have_ref "BC10 control scans (our pools, the stored templates)" "$TPL" || exit 0 ;;
  *) echo "unknown variant $VARIANT" >&2; exit 2 ;;
esac
OUT="$W/scan_$VARIANT"; mkdir -p "$OUT"
MASSES="${*:-$(tr '\n' ' ' < "$POINTS")}"
export ALP W OUT TPL VARIANT SLOT L

one_mass() {   # one HTCondor job of the original, run locally
  set -uo pipefail
  local M=$1 D="$OUT/bc10_$1" T="$OUT/bc10_$1.tmp" t0
  [ -s "$D/sensitivity.csv" ] && { echo "skip $M (done)"; return 0; }
  rm -rf "$T"; mkdir -p "$T/out" "$T/geo" "$T/work"
  t0=$(date +%s)
  GRENDEL_WORK_DIR="$T/work" MPLCONFIGDIR="$W/.mpl" \
  $SLOT -n 1 --label "alp:04_scan_${VARIANT}_$M" -- /usr/bin/time -f "%e s wall, %M kB max RSS" -o "$T/time.txt" \
    python -m grendel.models.alp.scan --mass "$M" --workers 1 --thresholds 3 10 \
      --vectors-dir "$W/llp_4vectors" --templates-dir "$TPL" --geometry-dir "$T/geo" --out "$T/out" \
      > "$T/scan.log" 2>&1 || { echo "FAIL $M (see $T/scan.log)"; return 1; }
  cp "$T"/out/*.csv "$T"/out/*.json "$T/" 2>/dev/null
  [ -s "$T/sensitivity.csv" ] || { echo "FAIL $M (no row)"; return 1; }
  rm -rf "$T/out" "$T/work"             # keep geo/ (the ray-cast cache, ~0.3 MB): validate_geometry.py reads it
  mv "$T" "$D"
  echo "done $M in $(( $(date +%s) - t0 )) s ($(cat "$D/time.txt")): $(tail -1 "$D/sensitivity.csv" | cut -d, -f1-6)"
}
export -f one_mass

echo "$(date '+%F %T') 04_scan $VARIANT start (J=$J): $(echo $MASSES | wc -w) masses" | tee -a "$L/04_scan_$VARIANT.log"
# shellcheck disable=SC2086
printf '%s\n' $MASSES | xargs -P "$J" -I{} bash -c 'one_mass {}' >> "$L/04_scan_$VARIANT.log" 2>&1 || true
n_done=$(ls -d "$OUT"/bc10_*/ 2>/dev/null | grep -vc '\.tmp/$' || true)
echo "$(date '+%F %T') 04_scan $VARIANT end: $n_done rows present, $(grep -c '^FAIL' "$L/04_scan_$VARIANT.log" || true) FAIL lines in the log" \
  | tee -a "$L/04_scan_$VARIANT.log"
