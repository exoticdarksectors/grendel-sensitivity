#!/usr/bin/env bash
# Step 01 -- BC4 four-vectors: FONLL bottom -> b-hadron (B+, B0, Bs, Lambda_b) -> X_s S, 86 masses.
#
# Published recipe: the v1 pool (curves/data/grendel_v2/README.md; scratch/bc4_nsig10_20260820/llp_4vectors, used
# unchanged by bc4_exhad_20260918, v2_baseline, v2.3 and the v2.5 scans): --n-pool 400000 --high-pt-tilt-scale 5
# (50 % nominal mixture), seed 42, the 86 masses of production.MASS_GRID in one RNG stream, the tracked FONLL grid,
# `particle` 0.26.2. This is grendel.models.scalar.scan's produce_missing (py/produce_bc4.py calls the same
# functions in the same order, writing each file atomically).
# Output: work/scalar/bc4/llp_4vectors/mS_<m>.csv (86 x 122 MB, 1.6e6 rows; 1.2e6 at 4.6-4.7 GeV).
# Validation: byte sizes against v2_3_20261016/lxplus/inputs_manifest.txt (the files themselves are HFS-compressed
# on the Mac drive and read as empty here); the contents are checked by the control scans of 07_scan_controls.sh.
# Without the reference data ($GRENDEL_REF) the size check is skipped.
# Resume-safe (complete files are kept). 1 slot, ~6 min, ~1 GB RAM.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
step_lock
V="$W/bc4/llp_4vectors"
n=$(find "$V" -maxdepth 1 -name 'mS_*.csv' 2>/dev/null | wc -l)
if [ "$n" -lt 86 ]; then   # (all present: nothing to produce, no slot needed)
  # shellcheck disable=SC2086
  PYTHONPATH="$PARTICLE_OVERLAY" $SLOT -n 1 --label scalar:01_bc4_vectors -- \
    python "$S/py/produce_bc4.py" "$V" "$BC4_N_POOL" "$PROD_SEED" "$BC4_TILT" $BC4_MASSES >> "$L/01_bc4_vectors.log" 2>&1
fi
mkdir -p "$W/validation"
if ref_ok "BC4 four-vector sizes vs the published pool" "$REF_BC4_MANIFEST"; then
  python "$S/py/validate_vectors.py" sizes "$V" "$REF_BC4_MANIFEST" bc4/llp_4vectors/ "$W/validation/bc4_vectors_sizes.tsv" \
    | tee "$W/validation/bc4_vectors.txt"
fi
step_ok 01_bc4_vectors
