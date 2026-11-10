#!/usr/bin/env bash
# Step 04 -- BC5 Higgs production, stage 3: the Higgs (pT, y) histogram per mode (normalised to YR4 by the
# production), as grendel/models/scalar/data/higgs_pt_y_14TeV.csv.gz.
#
# Published recipe: python -m grendel.models.scalar.higgs_pool (2026-10-14 17:48, from asbuilt_20261014/bc5/madgraph;
# the header of the tracked file records it). Here written to work/ instead of the package (src/ stays untouched):
#   (a) work/scalar/bc5/higgs_hist/higgs_pt_y_14TeV.csv.gz            from our MG5 + shower (02, 03a)
#   (b) work/scalar/bc5/higgs_hist/higgs_pt_y_14TeV_published_lhe.csv.gz  from the stored LHE + our shower (03b)
# Validation (validation/higgs_hist.txt): (b) against the tracked file cell by cell (the shower alone), (a) against it
# with a chi^2 per mode (MG5 statistics). (b) exists only where 03 had the reference data ($GRENDEL_REF); otherwise
# it is skipped. Seconds.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
step_lock
H="$W/bc5/higgs_hist"
TRACKED="$SRC/grendel/models/scalar/data/higgs_pt_y_14TeV.csv.gz"
mkdir -p "$H" "$W/validation"
[ -s "$H/higgs_pt_y_14TeV.csv.gz" ] || \
  python -m grendel.models.scalar.higgs_pool --work-dir "$W/bc5/madgraph" --out "$H/higgs_pt_y_14TeV.csv.gz" \
    >> "$L/04_bc5_higgs_hist.log" 2>&1
RS="$W/bc5/shower_of_published_lhe"
HAVE_B=0
if [ -s "$H/higgs_pt_y_14TeV_published_lhe.csv.gz" ]; then
  HAVE_B=1
elif ref_ok "(b) histogram of the stored LHE files showered here (03b)" "$RS/higgs_shower_manifest.json"; then
  python -m grendel.models.scalar.higgs_pool --work-dir "$RS" \
    --out "$H/higgs_pt_y_14TeV_published_lhe.csv.gz" >> "$L/04_bc5_higgs_hist.log" 2>&1
  HAVE_B=1
fi
{
  if [ "$HAVE_B" -eq 1 ]; then
    echo "== (b) published LHE, showered here, vs the tracked histogram"
    python "$S/py/compare_higgs_hist.py" "$H/higgs_pt_y_14TeV_published_lhe.csv.gz" "$TRACKED"
  fi
  echo "== (a) MG5 + shower regenerated here vs the tracked histogram"
  python "$S/py/compare_higgs_hist.py" "$H/higgs_pt_y_14TeV.csv.gz" "$TRACKED"
} | tee "$W/validation/higgs_hist.txt"
step_ok 04_bc5_higgs_hist
