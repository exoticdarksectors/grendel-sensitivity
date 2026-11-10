#!/usr/bin/env bash
# Step 02 -- BC5 Higgs production, stage 1: MadGraph5_aMC@NLO 3.6.6 LO hard process of ggF (heft), VBF, WH, ZH and
# ttH at 14 TeV, NNPDF40_nlo_as_01180, 200 000 unweighted events per mode (the Higgs undecayed).
#
# Published recipe: grendel-sensitivity/grendel_work/bc5/madgraph (2026-09-25, macOS): higgs_production
# --nevents 200000 with `generate_events ... --multicore --nb_core=4` (generate_events_run_200000.log), the run card
# of cards/run_card_template.dat (iseed 0: MG5 chose seed offset 21 in every fresh process folder).
#   GRENDEL_BC5_MG5_WORK_DIR=work/scalar/bc5/madgraph python -m grendel.models.scalar.higgs_production \
#       --nevents 200000 --nb-core 4
# MG5's own work dir is ours (GRENDEL_BC5_MG5_WORK_DIR); the global MG5 configuration is not touched.
# Output: work/scalar/bc5/madgraph/higgs_<mode>/Events/run_200000/unweighted_events.lhe.gz + higgs_lhe_manifest.json.
# Validation: LO cross sections and random seeds against the stored banners (validation/mg5.txt); the kinematics
# are compared after the shower (04); skipped without the reference data ($GRENDEL_REF). Resume-safe: modes already
# in the manifest with their LHE are skipped.
# 4 slots, ~1-2 h.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
step_lock
K=${SCALAR_MG5_SLOTS:-4}
export GRENDEL_BC5_MG5_WORK_DIR="$W/bc5/madgraph"
mkdir -p "$GRENDEL_BC5_MG5_WORK_DIR" "$W/validation"
todo=$(python - <<'EOF'
import json, os
from pathlib import Path
w = Path(os.environ["GRENDEL_BC5_MG5_WORK_DIR"]); m = w / "higgs_lhe_manifest.json"
done = json.loads(m.read_text()) if m.exists() else {}
print(" ".join(x for x in ("ggf", "vbf", "wh", "zh", "tth") if not (x in done and Path(done[x]["lhe"]).exists())))
EOF
)
if [ -n "$todo" ]; then
  echo "$(date '+%F %T') 02_bc5_higgs_mg5: modes $todo" | tee -a "$L/02_bc5_higgs_mg5.log"
  # shellcheck disable=SC2086
  $SLOT -n "$K" --label scalar:02_mg5 -- \
    python -m grendel.models.scalar.higgs_production --modes $todo --nevents "$MG5_NEVENTS" --nb-core "$K" \
    >> "$L/02_bc5_higgs_mg5.log" 2>&1
fi
# validation (skipped when its report is newer than the LHE manifest, or without the stored run)
if ! ref_ok "MG5 cross sections, seeds and LHE vs the stored 2026-09-25 run" "$REF_MG5"; then
  :
elif [ "$W/validation/mg5.txt" -nt "$GRENDEL_BC5_MG5_WORK_DIR/higgs_lhe_manifest.json" ] && grep -q "tth LHE" "$W/validation/mg5.txt"; then
  cat "$W/validation/mg5.txt"
else
  python - "$GRENDEL_BC5_MG5_WORK_DIR" "$REF_MG5" <<'EOF' | tee "$W/validation/mg5.txt"
import gzip, json, re, sys
from pathlib import Path
ours, ref = Path(sys.argv[1]), Path(sys.argv[2])
def banner(d, mode):
    b = next((d / f"higgs_{mode}" / "Events" / "run_200000").glob("*banner.txt")).read_text()
    xs = re.search(r"Integrated weight \(pb\)\s*:\s*([0-9.eE+-]+)", b).group(1)
    seed = re.search(r"(\d+)\s*=\s*iseed", b).group(1)
    return float(xs), int(seed)
print("mode  sigma_LO ours [pb]  published  ratio  seed ours/published")
for mode in ("ggf", "vbf", "wh", "zh", "tth"):
    (xo, so), (xr, sr) = banner(ours, mode), banner(ref, mode)
    print(f"{mode:4s}  {xo:12.6g}  {xr:12.6g}  {xo / xr:.5f}  {so}/{sr}")
EOF
  for mode in ggf vbf wh zh tth; do
    echo -n "$mode LHE vs published: " | tee -a "$W/validation/mg5.txt"
    python "$S/py/compare_lhe.py" "$GRENDEL_BC5_MG5_WORK_DIR/higgs_$mode/Events/run_200000/unweighted_events.lhe.gz" \
      "$REF_MG5/higgs_$mode/Events/run_200000/unweighted_events.lhe.gz" | tee -a "$W/validation/mg5.txt"
  done
fi
step_ok 02_bc5_higgs_mg5
