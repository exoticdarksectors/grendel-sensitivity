#!/usr/bin/env bash
# 23_check_scalar_1809.sh -- (d) Mathematica-side cross-check of the scalar decay tables behind Figs. 4 and 5.
#   source env/activate.sh && bash pipelines/senscalc/23_check_scalar_1809.sh
# grendel/models/scalar/data/eventcalc_1809/{ctau,branching_ratios}.json.gz are copied from exHad (scalar-1809), which took
# them from EventCalc/SensCalc (PROVENANCE.md there). The SensCalc originals are MX dumps,
# phenomenology/Scalar/decay widths/{ctauS,BrRatioScalar}.mx (identical in v1.4 and eb2c456). This decodes ctauS.mx on
# Linux and compares the lifetime with the tracked table; BrRatioScalar.mx is decoded to its structure only.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
OUT="$SCW/scalar_1809"; mkdir -p "$OUT"
for f in ctauS.mx BrRatioScalar.mx; do
  a=$(sha256sum < "$REPRO/third_party/SensCalc/phenomenology/Scalar/decay widths/$f" | cut -c1-16)
  b=$(sha256sum < "$REPRO/third_party/SensCalc-v1.4/phenomenology/Scalar/decay widths/$f" | cut -c1-16)
  echo "$f eb2c456 $a v1.4 $b $([ "$a" = "$b" ] && echo same || echo DIFFERENT)"
done | tee "$SCW/validation/scalar_1809.txt"
if [ ! -f "$MARK/decode_scalar_1809" ]; then
  require_kernel
  with_slots 1 23_decode wolframscript -file "$WL/decode_scalar_1809.wls" "$REPRO/third_party/SensCalc" "$OUT" < /dev/null 2>&1 | tee "$SCL/23_decode_scalar_1809.log"
  [ -f "$OUT/ctauS_1.csv" ] || { echo "decode failed"; exit 1; }
  touch "$MARK/decode_scalar_1809"
fi
python3 -I "$PY/compare_scalar_ctau.py" "$OUT/ctauS_1.csv" \
  "$REPRO/src/grendel-sensitivity/grendel/models/scalar/data/eventcalc_1809/ctau.json.gz" | tee -a "$SCW/validation/scalar_1809.txt"
