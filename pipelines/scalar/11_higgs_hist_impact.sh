#!/usr/bin/env bash
# Step 11 -- how much the BC5 rows depend on which Higgs histogram is used: the tracked one (the published input,
# used by 05) or the one rebuilt here from scratch (02-04: MG5 + Pythia on Linux, statistically equivalent but not
# identical, see README).
#
# Same RNG stream for both (seed 42, the scan_bc5 defaults, masses 1 15 50 in their own stream, so the rows are not
# the production rows), so the only difference is the histogram; each set is scanned with the published settings and
# our templates. 1 GeV: the h -> SS rows set the upper edge; 15 and 50 GeV: h -> SS is the only channel.
# Output: work/scalar/bc5/hist_impact/{tracked,rebuilt}/ and validation/higgs_hist_impact.txt.
# ~10 min of CPU (2 x 3 small productions, 6 scans), 1 slot per job. Resume-safe.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
step_lock
MASSES="1 15 50"
HI="$W/bc5/hist_impact"
for v in tracked rebuilt; do
  hist=tracked; [ $v = rebuilt ] && hist="$W/bc5/higgs_hist/higgs_pt_y_14TeV.csv.gz"
  # shellcheck disable=SC2086
  $SLOT -n 1 --label scalar:11_prod_$v -- \
    python "$S/py/produce_bc5.py" "$HI/$v/llp_4vectors" "$BC5_N_POOL" "$BC5_N_HIGGS" "$PROD_SEED" "$hist" $MASSES \
    >> "$L/11_higgs_hist_impact.log" 2>&1
  for m in $MASSES; do
    O="$HI/$v/out/$m"
    [ -s "$O/time.txt" ] && continue
    mkdir -p "$O"
    $SLOT -n 1 --label "scalar:11_scan_${v}_$m" -- bash "$S/scan_job.sh" "$O" \
      python -m grendel.models.scalar.scan_bc5 --workers 1 --decay-samples 100 --thresholds 3 10 \
        --vectors-dir "$HI/$v/llp_4vectors" --templates-dir "$W/templates_exhad" --geometry-dir "$HI/$v/geo" \
        --out "$O" --mass "$m" >> "$L/11_higgs_hist_impact.log" 2>&1 &
  done
done
wait
python - "$HI" <<'EOF' | tee "$W/validation/higgs_hist_impact.txt"
import sys, warnings
from pathlib import Path
import pandas as pd
warnings.filterwarnings("ignore")
hi = Path(sys.argv[1])
print("BC5 rows with the rebuilt Higgs histogram / with the tracked one (same RNG stream, same templates):")
for m in ("1", "15", "50"):
    a = pd.read_csv(hi / "rebuilt/out" / m / "sensitivity.csv").iloc[0]
    b = pd.read_csv(hi / "tracked/out" / m / "sensitivity.csv").iloc[0]
    print(f"  {m:>3} GeV: u2_min x{a.u2_min / b.u2_min:.4f}, u2_max x{a.u2_max / b.u2_max:.4f}, "
          f"peak N x{a.peak_N / b.peak_N:.4f}; N>=10 edges x{a.u2_min_N10 / b.u2_min_N10:.4f}, "
          f"x{a.u2_max_N10 / b.u2_max_N10:.4f}")
EOF
echo "11_higgs_hist_impact OK"
