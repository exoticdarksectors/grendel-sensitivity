#!/usr/bin/env bash
# The whole BC10 chain at the published statistics, resume-safe (every step skips what is done).
#
#   source env/activate.sh
#   setsid nohup bash pipelines/alp/run_all.sh > logs/alp/run_all.out 2>&1 < /dev/null &
#   cat logs/alp/run_all.status            # one line per finished step
#
# Knobs (environment): ALP_K = parallel one-slot template jobs (default 4), ALP_J = main scan jobs in flight
# (default 5), ALP_JC = control scan jobs in flight (default 2; they run alongside the templates). Every CPU-heavy
# command goes through pipelines/common/slot_run.py with one slot per process.
# Measured on this laptop (Ryzen 5 PRO 7530U): 02 ~7 min, 03 ~50 s per mass per core, 04 ~90-100 s per mass.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
K="${ALP_K:-4}"; J="${ALP_J:-5}"; JC="${ALP_JC:-2}"
ST="$L/run_all.status"
step() { local name=$1; shift; local t0; t0=$(date +%s)
  if "$@"; then echo "$(date '+%F %T') OK    $name ($(( $(date +%s) - t0 )) s)" >> "$ST"
  else echo "$(date '+%F %T') FAIL  $name ($(( $(date +%s) - t0 )) s)" >> "$ST"; return 1; fi; }

echo "$(date '+%F %T') START run_all (K=$K, J=$J, JC=$JC, pid $$)" >> "$ST"
step 01_cbs bash "$ALP/01_cbs.sh" > /dev/null
step 02_production bash "$ALP/02_production.sh"
# byte identity of the pools: expected to FAIL on this platform (README, "02, four-vectors"), so not fatal; the
# validate_* scripts print SKIP and exit 0 without the reference tree ($GRENDEL_REF)
step "validate_pools (byte identity, expected to differ)" \
  python "$ALP/validate_pools.py" "$W/llp_4vectors" "$W/validation/pools.txt" || true

# control branch: our pools with the stored macOS templates of the published run (needs no step-03 output; skipped
# without the reference tree, $GRENDEL_REF)
control=""
if have_ref "BC10 control branch (04 control, 05 control, validate_geometry)" "$CAMP/templates_exhad"; then
  (
    step "04_scan control" bash "$ALP/04_scan.sh" control "$JC" &&
    step "05_collect control" python "$ALP/05_collect.py" control > /dev/null &&
    step validate_geometry python "$ALP/validate_geometry.py" "$W/scan_control" "$W/validation/geometry.txt"
  ) &
  control=$!
else
  echo "$(date '+%F %T') SKIP  control branch (no reference tree)" >> "$ST"
fi

# main branch: our templates
step 03_templates bash "$ALP/03_templates.sh" "$K"
step validate_templates python -W ignore "$ALP/validate_templates.py" "$W/templates_exhad" "$W/validation/templates.txt"
step "04_scan main" bash "$ALP/04_scan.sh" main "$J"
step "05_collect main" python "$ALP/05_collect.py" main > /dev/null
step 06_plot bash "$ALP/06_plot.sh" main > /dev/null
[ -z "$control" ] || wait "$control" || true
echo "$(date '+%F %T') END run_all" >> "$ST"
