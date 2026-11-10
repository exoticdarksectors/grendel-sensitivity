#!/usr/bin/env bash
# Step 06 -- exHad scalar-1809 rest-frame decay templates for every BC4 and BC5 mass (112; BC4's 86 are a subset).
#
# Published recipe: v2_baseline_20260929/templates/scalar (exHad 7bb2308, Pythia 8.317, 20 000 per mass, base seed
# 1234, sub-requests of 2000 with up to 3 retries each, chunk 512): the 86 BC4 masses from bc4_exhad_20260918
# (grendel-mono's generator, which the committed one reproduces bit for bit: v2_baseline_20260929/verify_port) and
# the 26 masses 4.8-62 GeV from logs/templates_bc5.log:
#   $EXHAD_PY -m grendel.models.scalar.templates_exhad --mass <MASS_GRID_BC5> --n-templates 20000 --seed 1234 \
#       --workers 6 --skip-existing --out templates/scalar
# exHad is the pinned e41b07b: same scalar data as 7bb2308 and the same scalar golden digests upstream (README).
# Every point and sub-request has its own seed, so the result does not depend on --workers.
# Output: work/scalar/templates_exhad/templates_<m>.npz (112, ~5 MB each) + MANIFEST.json; validation report
# work/scalar/validation/templates.tsv (skipped without the reference data, $GRENDEL_REF; 09 and 10 then cannot
# single out the masses whose decays differ from the published bundles). Resume-safe (--skip-existing). 4 slots, ~5 min.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
step_lock
K=${SCALAR_TPL_SLOTS:-4}
T="$W/templates_exhad"
MASSES=$(python -c "from grendel.models.scalar.production_bc5 import MASS_GRID_BC5 as G; print(' '.join(f'{m:g}' for m in G))")
n=$(find "$T" -maxdepth 1 -name "templates_*.npz" 2>/dev/null | wc -l)
if [ "$n" -lt 112 ]; then
  echo "$(date '+%F %T') 06_templates: generating ($n/112 present)" | tee -a "$L/06_templates.log"
  # shellcheck disable=SC2086
  PYTHONPATH="$SRC" $SLOT -n "$K" --label scalar:06_templates -- \
    "$GRENDEL_EXHAD_PYTHON" -m grendel.models.scalar.templates_exhad --model scalar-1809 --mass $MASSES \
      --n-templates "$N_TEMPLATES" --seed "$TPL_SEED" --workers "$K" --skip-existing --out "$T" \
      >> "$L/06_templates.log" 2>&1
fi
n=$(find "$T" -maxdepth 1 -name "templates_*.npz" | wc -l)
[ "$n" -eq 112 ] || { echo "06_templates FAILED: $n/112"; exit 1; }
mkdir -p "$W/validation"
# (validation skipped when its report is newer than every bundle, or without the published bundles)
if ref_ok "exHad templates vs the published bundles" "$REF_TEMPLATES"; then
  if [ -n "$(find "$T" -name 'templates_*.npz' -newer "$W/validation/templates.txt" 2>/dev/null | head -n 1)$( [ -s "$W/validation/templates.txt" ] || echo missing)" ]; then
    $SLOT -n 1 --label scalar:06_validate -- \
      python "$S/py/validate_templates.py" "$T" "$REF_TEMPLATES" "$W/validation/templates.tsv" | tee "$W/validation/templates.txt"
  else
    cat "$W/validation/templates.txt"
  fi
fi
step_ok 06_templates
