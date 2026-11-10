#!/bin/bash
# HNL step 07: exHad rest-frame decay templates for the 241 masses per flavour of the table: the recipe of the
# published sets (v2_baseline_20260929/templates/hnl = hnl_exhad_20260917, and hnl_fine_condor_20261014/templates.sh):
#   grendel.models.hnl.templates.exhad --n-templates 20000 --seed 1234 (chunk 512), exHad lifetime, central variation.
# exHad is the pin e41b07b (the published sets used 7bb2308; e41b07b changes only the ALP model: every HNL golden
# digest of exHad's test suite is unchanged between the two). Not done: --reference-templates (the FairShip ctau
# stored for diagnostics as ctau_m_u2eq1_reference; the scan does not read it).
# Six single-slot jobs (flavour x coarse/fine), each into its own directory work/hnl/templates_gen/<fl>_<grid>
# (own MANIFEST.json), linked into work/hnl/templates/<fl>/templates_<label>.npz. Resume-safe (--skip-existing).
set -euo pipefail
source "$(dirname "$0")/../../env/activate.sh"
T=$REPRO/work/hnl/templates_gen; mkdir -p $T $REPRO/work/hnl/templates
L=$REPRO/logs/hnl/07; mkdir -p $L
export PYTHONPATH=$REPRO/src/grendel-sensitivity PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 MPLCONFIGDIR=$T/.mpl
masses() { python -c "import sys; sys.path.insert(0, '$REPRO/pipelines/hnl/lib'); import grids; print(*grids.$1)"; }
gen() {  # gen <flavor> <COARSE|FINE>   (flock: one generator per directory, also across invocations)
  local out=$T/${1}_$2
  flock $T/.lock_${1}_$2 python $REPRO/pipelines/common/slot_run.py -n 1 --label hnl:07_templates_${1}_$2 -- \
    $GRENDEL_EXHAD_PYTHON -m grendel.models.hnl.templates.exhad --flavor $1 --mass $(masses $2) \
      --n-templates 20000 --seed 1234 --chunk-size 512 --workers 1 --skip-existing \
      --exhad-root $REPRO/third_party/exHad --out $out >> $L/${1}_$2.log 2>&1
}
# work/hnl/templates/<fl>/templates_<label>.npz -> the generator's file (linked before generation: a scan waits
# until the target exists)
python - <<'PY'
import os, sys
from pathlib import Path
sys.path.insert(0, os.environ["REPRO"] + "/pipelines/hnl/lib")
import grids
W = Path(os.environ["REPRO"]) / "work/hnl"
for fl in grids.FLAVORS:
    (W / "templates" / fl).mkdir(parents=True, exist_ok=True)
    for g, ms in (("COARSE", grids.COARSE), ("FINE", grids.FINE)):
        for m in ms:
            name = f"templates_{grids.label(m)}.npz"
            link = W / "templates" / fl / name
            if not link.is_symlink():
                link.symlink_to(W / "templates_gen" / f"{fl}_{g}" / fl / name)
PY
# a generator killed while writing leaves a truncated bundle that --skip-existing would keep: drop unreadable ones
$GRENDEL_EXHAD_PYTHON - $T <<'PY'
import sys, zipfile
from pathlib import Path
import numpy as np
for p in Path(sys.argv[1]).glob("*/*/templates_*.npz"):
    try:
        with np.load(p) as d:
            d["daughter_counts"], d["energy"], d["ctau_m_u2eq1"]
    except (OSError, ValueError, KeyError, EOFError, zipfile.BadZipFile) as err:
        print(f"removing unreadable {p}: {err}")
        p.unlink()
PY
pids=()
for fl in Ue Umu Utau; do for g in COARSE FINE; do gen $fl $g & pids+=($!); done; done
rc=0; for p in "${pids[@]}"; do wait $p || rc=1; done
n=$(ls $REPRO/work/hnl/templates/*/templates_*.npz | wc -l)
echo "templates: $n of 723"
[ $rc -eq 0 ] && [ $n -eq 723 ]
