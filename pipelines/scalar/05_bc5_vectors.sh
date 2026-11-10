#!/usr/bin/env bash
# Step 05 -- BC5 four-vectors: b -> X_s S (mixing), h -> S S, b -> X_s S S and B_s -> S S, 112 masses.
#
# Published recipe: asbuilt_20261014/chain.sh (code 461b7e5, = 8802f50 for the production), the scan_bc5 defaults:
# --n-pool 200000 --n-higgs 200000 --seed 42, BR(h->SS) 0.01, quartic B on, the tracked FONLL grid, the tracked
# SensCalc quartic tables and the tracked Higgs (pT, y) histogram, one RNG stream (three generators) over the
# 112 masses in the interleaved order of BC5_MASSES. py/produce_bc5.py calls scan_bc5's produce_missing functions
# in the same order and writes each file atomically.
# HIGGS_HIST (default `tracked`, the published input) can point at the histogram rebuilt by 02-04
# (work/scalar/bc5/higgs_hist/higgs_pt_y_14TeV.csv.gz); see the README for why the default is the tracked one.
# Output: work/scalar/bc5/llp_4vectors/mS_<m>.csv (112 x 181 MB). Validation against the stored vectors
# asbuilt_20261014/bc5/llp_4vectors (SHA-256, else the differing lines and their largest relative difference):
# work/scalar/validation/bc5_vectors.tsv; skipped without the reference data ($GRENDEL_REF). Resume-safe. 1 slot,
# ~15 min, ~1.5 GB RAM.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
step_lock
HIGGS_HIST=${HIGGS_HIST:-tracked}
V="$W/bc5/llp_4vectors"
n=$(find "$V" -maxdepth 1 -name 'mS_*.csv' 2>/dev/null | wc -l)
if [ "$n" -lt 112 ]; then   # (all present: nothing to produce, no slot needed)
  # shellcheck disable=SC2086
  $SLOT -n 1 --label scalar:05_bc5_vectors -- \
    python "$S/py/produce_bc5.py" "$V" "$BC5_N_POOL" "$BC5_N_HIGGS" "$PROD_SEED" "$HIGGS_HIST" $BC5_MASSES \
    >> "$L/05_bc5_vectors.log" 2>&1
fi
mkdir -p "$W/validation"
# (validation skipped when its report is newer than every four-vector file, or without the stored vectors)
if [ "$HIGGS_HIST" = tracked ] && ref_ok "BC5 four-vectors vs the stored pool" "$REF_BC5_VECTORS" \
   && [ -n "$(find "$V" -name 'mS_*.csv' -newer "$W/validation/bc5_vectors.txt" 2>/dev/null | head -n 1)$( [ -s "$W/validation/bc5_vectors.txt" ] || echo missing)" ]; then
  $SLOT -n 1 --label scalar:05_validate -- \
    python "$S/py/validate_vectors.py" bytes "$V" "$REF_BC5_VECTORS" "$W/validation/bc5_vectors.tsv" \
    | tee "$W/validation/bc5_vectors.txt"
fi
step_ok 05_bc5_vectors
