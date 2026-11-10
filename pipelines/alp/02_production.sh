#!/usr/bin/env bash
# Step 02 -- BC10 ALP four-vectors (B -> K^(i) a), the v2.2 recipe at full statistics.
#
#   (a) v2.2 pool  = bc10_gkoz_cbs_20261005 (curves/data/grendel_v2_2/README.md): the default BC10 grid
#       (ALP_MASS_GRID, 145 masses, eta/eta' skipped) in ONE RNG stream, seed 42:
#         python -m grendel.models.alp.production --n-pool 1200000 --high-pt-tilt-scale 5 --out-dir ...
#   (b) tip pool   = bc10_v24_20261019/lxplus/job_prod.sh: its own stream over
#         --mass 3.4 3.42 3.44 3.46 3.48 3.5 3.52 3.54 3.56 3.58 (same options); the 8 tip masses are used.
#   (c) work/alp/llp_4vectors: links to (a) for the 145 grid masses and to (b) for the 8 tip masses, as the
#       v2.4 inputs/llp_4vectors on AFS (links to the v2.2/v2.3 vectors + 8 tip files).
#
# The production draws every mass from one RNG stream, so a mass is reproduced only inside the full run with the
# same mass list; --resume replays the stream and keeps finished files (atomic writes), hence resume-safe.
# Outputs: work/alp/vectors_v22/mA_*.csv (143 x ~92 MB), work/alp/vectors_tip/, work/alp/llp_4vectors/.
# Two slots (the two runs in parallel, one core each); ~1 GB RAM each.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
OPTS=(--n-pool "$N_POOL" --high-pt-tilt-scale "$TILT" --seed "$PROD_SEED" --resume)

run_v22() {
  $SLOT -n 1 --label alp:02_prod_v22 -- \
    python -m grendel.models.alp.production "${OPTS[@]}" --out-dir "$W/vectors_v22" >> "$L/02_prod_v22.log" 2>&1
}
run_tip() {
  # shellcheck disable=SC2086
  $SLOT -n 1 --label alp:02_prod_tip -- \
    python -m grendel.models.alp.production --mass $TIP_RUN_MASSES "${OPTS[@]}" --out-dir "$W/vectors_tip" \
    >> "$L/02_prod_tip.log" 2>&1
}
done_v22() { [ "$(grep -c '' <(ls "$W/vectors_v22"/mA_*.csv 2>/dev/null))" -eq 143 ]; }   # 145 grid masses, eta and eta' skipped
done_tip() { for m in $TIP_RUN_MASSES; do [ -f "$W/vectors_tip/mA_$(printf '%.3f' "$m" | tr . p).csv" ] || return 1; done; }

echo "$(date '+%F %T') 02_production start" | tee -a "$L/02_prod.log"
pids=()
done_v22 || { run_v22 & pids+=($!); }
done_tip || { run_tip & pids+=($!); }
rc=0; for p in "${pids[@]}"; do wait "$p" || rc=1; done
[ $rc -eq 0 ] && done_v22 && done_tip || { echo "02_production FAILED (see $L/02_prod_*.log)"; exit 1; }

# (c) the scan input folder
mkdir -p "$W/llp_4vectors"
for f in "$W/vectors_v22"/mA_*.csv; do ln -sfn "../vectors_v22/$(basename "$f")" "$W/llp_4vectors/$(basename "$f")"; done
for m in $TIP_MASSES; do
  b="mA_$(printf '%.3f' "$m" | tr . p).csv"; ln -sfn "../vectors_tip/$b" "$W/llp_4vectors/$b"
done
echo "$(date '+%F %T') 02_production done: $(ls "$W/llp_4vectors" | wc -l) four-vector files" | tee -a "$L/02_prod.log"
