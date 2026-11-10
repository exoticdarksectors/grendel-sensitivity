# lib.sh -- helpers sourced by the numbered scripts of pipelines/senscalc (not run directly).
: "${REPRO:?source env/activate.sh first}"
set -euo pipefail

SCP="$REPRO/pipelines/senscalc"          # this pipeline
WL="$SCP/wl"                             # Wolfram drivers (campaign copies + senscalc_env.wl)
PY="$SCP/py"
SCW="$REPRO/work/senscalc"               # large outputs
SCL="$REPRO/logs/senscalc"
MARK="$SCW/markers"
EB=eb2c456e8c9b35e77d212d18813306da195a13cb
# Reference material (read only, optional: $REF of pipelines/common/ref.sh; without it the comparisons are skipped):
#   the BC10 SensCalc campaign the published FASER2-FPF curve comes from (its SensCalc tree with every output, logs)
#   and M. Ovchynnikov's SHiP JSON (private communication 2026-10-16, not public), the source of the tracked
#   curves/data/bc10/projections/ship_senscalc.dat.
source "$REPRO/pipelines/common/ref.sh"
CAMPAIGN=$REF/grendel/campaigns/bc10/faser2_senscalc_20261025
CAMPAIGN_TREE="$CAMPAIGN/src/SensCalc@eb2c456"
MAKSYM_JSON=$REF/grendel/campaigns/bc10/maksym_alp_fix_check_20261018/ship_json/Sensitivity_ALP-fermion_at_SHiP-ECN3_Nev=2.3_Npot=6.e20_2026-10-16.json
NK="${SENSCALC_SUBKERNELS:-4}"          # parallel subkernels per SensCalc run (= slots taken); the campaign used 4
mkdir -p "$SCW" "$SCL" "$MARK"

# ref_or_skip REPORT LABEL PATH... -- have_ref (pipelines/common/ref.sh) that also writes its SKIP line into REPORT
ref_or_skip() {
  local rep=$1 m; shift
  if m=$(have_ref "$@" 2>&1); then return 0; fi
  echo "$m" | tee -a "$rep" >&2
  return 1
}

# the kernel must be activated (wolframscript -activate, with the user's Wolfram ID)
require_kernel() {
  local out
  out=$(timeout 120 wolframscript -code 'Print["OK ", $VersionNumber]' < /dev/null 2>&1 || true)
  case "$out" in
    *"OK "*) ;;
    *) echo "ERROR: the Wolfram Engine is not usable: $out"
       echo "       activate it once with:  source env/activate.sh && wolframscript -activate"; exit 3 ;;
  esac
}

# wait until enough memory is available (SensCalc's main kernel peaked at 4.3 GB on macOS)
wait_mem() {
  local need_gb=${1:-5}
  while [ "$(free -g | awk '/^Mem:/{print $7}')" -lt "$need_gb" ]; do
    echo "$(date '+%F %T') waiting for ${need_gb} GB of available memory"; sleep 60
  done
}

# Slots of the shared pool: every step goes through pipelines/common/slot_run.py -n K (first come, first served since
# 2026-11-06 21:53), K = cores the step uses (SensCalc runs: main kernel + $NK subkernels -> NK slots; light steps: 1).
#   SENSCALC_SPARE_CORE=1: a 1-core step under ~5 min (engine check, MX decoding) runs directly on the core the pool leaves
#   free (12 cores, 11 slots). Nothing long runs outside the pool.
with_slots() {
  local k="$1" label="$2"; shift 2
  if [ "$k" = 1 ] && [ "${SENSCALC_SPARE_CORE:-0}" = 1 ]; then "$@"; return; fi
  python3 "$REPRO/pipelines/common/slot_run.py" -n "$k" --label "senscalc:$label" -- "$@"
}

# run_wl <marker> <log> <slots> <wolfram driver> args...   -- one SensCalc notebook run, resume-safe
run_wl() {
  local marker="$MARK/$1" log="$2" k="$3"; shift 3
  if [ -f "$marker" ]; then echo "skip (done): $(basename "$marker")"; return 0; fi
  wait_mem 5
  echo "$(date '+%F %T') start $(basename "$marker") -> $log"
  local t0=$SECONDS
  with_slots "$k" "$(basename "$marker")" wolframscript -file "$@" > "$log" 2>&1 < /dev/null || true
  local dt=$((SECONDS - t0))
  if grep -q 'EMU| ABORTED\|UNANSWERED\|No valid password' "$log" || ! grep -q '^done ' "$log"; then
    echo "FAILED after ${dt} s: $(basename "$marker") (see $log)"; return 1
  fi
  local run; run=$(grep -o 'runMaster finished, total [0-9]* s' "$log" | grep -o '[0-9]* s' || echo "?")
  echo "${run} run, ${dt} s with queue" > "$marker"
  echo "$(date '+%F %T') done $(basename "$marker"): ${run} run, ${dt} s with queue"
}

# make_tree <dir>: copy of SensCalc eb2c456 (git archive, as the campaign's GitHub tarball) + the case aliases
make_tree() {
  local dir="$1"
  if [ -f "$dir/.tree_ok" ]; then python3 -I "$PY/case_links.py" "$dir" --apply > /dev/null; return 0; fi  # idempotent
  rm -rf "$dir"; mkdir -p "$dir"
  git -C "$REPRO/third_party/SensCalc" archive "$EB" | tar -x -C "$dir"
  # case-sensitive ext4 vs case-insensitive APFS: alias "phenomenology/<model>/production probabilities" (py/case_links.py)
  python3 -I "$PY/case_links.py" "$dir" --apply | tail -1
  echo "$EB" > "$dir/.tree_ok"
  echo "SensCalc tree ready: $dir ($(du -sh "$dir" | cut -f1))"
}
