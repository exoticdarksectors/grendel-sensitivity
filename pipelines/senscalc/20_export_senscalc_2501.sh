#!/usr/bin/env bash
# 20_export_senscalc_2501.sh -- re-export the tracked BC10 decay tables grendel/models/alp/data/senscalc_2501
# from the SensCalc v1.4 MX dumps with the repository's own exporter (export_senscalc_2501.py -> .wls), unchanged.
#   source env/activate.sh && bash pipelines/senscalc/20_export_senscalc_2501.sh
# Recipe: grendel/models/alp/data/senscalc_2501/PROVENANCE.md ("Regeneration"), commit db79c0e (2026-10-13; the
# exporter is unchanged from db79c0e to 8802f50). Source: third_party/SensCalc-v1.4 (4f7210b), hashes checked by the tool.
# Output: work/senscalc/senscalc_2501/{widths_raw_senscalc.csv, widths_bnt.csv, widths_metadata.json,
#         branching_ratios.csv, decay_channels.json, matrix_elements.json, EXPORT_MANIFEST.json} and *.gz (py/gzip_zlib.py)
# Validation: py/compare_export.py against the tracked files (byte-identical, else numerical with rel. tol 1e-12).
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
OUT="$SCW/senscalc_2501"
if [ ! -f "$MARK/export_senscalc_2501" ]; then
  require_kernel
  rm -rf "$OUT"
  ( cd "$REPRO/src/grendel-sensitivity" && with_slots 1 20_export \
      python -m grendel.models.alp.tools.export_senscalc_2501 "$REPRO/third_party/SensCalc-v1.4" --output-dir "$OUT" ) \
      2>&1 | tee "$SCL/20_export_senscalc_2501.log"
  [ -f "$OUT/EXPORT_MANIFEST.json" ] || { echo "export failed"; exit 1; }
  touch "$MARK/export_senscalc_2501"
fi
# the tracked .gz were written by macOS `gzip -n -9` (zlib); GNU gzip gives other bytes for the same content
for f in widths_bnt.csv branching_ratios.csv decay_channels.json matrix_elements.json; do
  python3 -I "$PY/gzip_zlib.py" "$OUT/$f" "$OUT/$f.gz"
done
python3 -I "$PY/compare_export.py" "$OUT" "$REPRO/src/grendel-sensitivity/grendel/models/alp/data/senscalc_2501" \
  | tee "$SCW/validation/senscalc_2501.txt"
