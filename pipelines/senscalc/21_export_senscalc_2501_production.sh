#!/usr/bin/env bash
# 21_export_senscalc_2501_production.sh -- decode the SensCalc v1.4 LHC production inputs of BC10 with the repository's
# export_senscalc_2501_production.py (-> .wls), unchanged. Its default destination is $GRENDEL_WORK_DIR/bc10/...;
# here it writes to work/senscalc/senscalc_2501_production. Nothing tracked is derived from it (BC10 production audit).
#   source env/activate.sh && bash pipelines/senscalc/21_export_senscalc_2501_production.sh
# Validation: the tool's own hash/manifest checks; informational comparison with the only stored export (optional
# reference tree), grendel/signal-models/bc10-alp/llpatcolliders_BC10_PR18/alp_fermion/tmp/senscalc_2501_production
# (macOS, 2026-07-13, SensCalc v.1.3.3 0bca050, older exporter): the coefficient, fragmentation-probability, Drell-Yan
# and bremsstrahlung inputs are byte-identical between 0bca050 and v1.4, the light-meson file and the notebooks are not.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
OUT="$SCW/senscalc_2501_production"
OLD=$REF/grendel/signal-models/bc10-alp/llpatcolliders_BC10_PR18/alp_fermion/tmp/senscalc_2501_production
if [ ! -f "$MARK/export_senscalc_2501_production" ]; then
  require_kernel
  rm -rf "$OUT"
  ( cd "$REPRO/src/grendel-sensitivity" && with_slots 1 21_export_prod \
      python -m grendel.models.alp.tools.export_senscalc_2501_production "$REPRO/third_party/SensCalc-v1.4" --output-dir "$OUT" ) \
      2>&1 | tee "$SCL/21_export_senscalc_2501_production.log"
  [ -f "$OUT/PRODUCTION_EXPORT_MANIFEST.json" ] || { echo "export failed"; exit 1; }
  touch "$MARK/export_senscalc_2501_production"
fi
if have_ref "SensCalc v1.4 production export vs the stored July 2026 export" "$OLD"; then
  python3 -I "$PY/compare_export.py" "$OUT" "$OLD" | grep -v MANIFEST | tee "$SCW/validation/senscalc_2501_production_vs_july.txt" || true
fi
