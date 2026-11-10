#!/usr/bin/env bash
# 15_noise_faser2.sh -- run-to-run Monte Carlo spread of the FASER2-FPF curve (validation only).
#   source env/activate.sh && bash pipelines/senscalc/15_noise_faser2.sh
# SensCalc does not seed its samplers (no SeedRandom in notebook 2), so every run draws a new sample. Below ~1 GeV the
# FASER2 edges come from B -> K(*) a (100 % of the lower edge, 70-90 % of the upper edge at 0.3-0.5 GeV), and FASER2
# sees only the very forward tail of that sample (theta < 2.3 mrad: 0.4-0.8 % of the B-decay ALPs). This step redraws
# only that input: a copy of the production tree, notebook 2 for LHC/B (10^6 events, as production), then notebook 3
# for FASER2-FPF (N >= 2.3 and 3). Acceptances, the other spectra and the deterministic bremsstrahlung are kept.
# Output: work/senscalc/noise3/SensCalc@eb2c456/Sensitivity domains/ALP-fermion/FASER2-FPF/*.json,
#         work/senscalc/validation/noise_faser2.txt (lower/upper edges: noise run vs production run vs campaign;
#         the campaign comparisons need the optional reference tree)
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
require_kernel
P="$SCW/run/SensCalc@eb2c456"
T="$SCW/noise3/SensCalc@eb2c456"
[ -f "$MARK/nb3_FASER2-FPF" ] || { echo "ERROR: production nb3_FASER2-FPF has not run"; exit 1; }
if [ ! -f "$T/.tree_ok" ]; then
  mkdir -p "$(dirname "$T")"; cp -a "$P" "$T.tmp" && mv "$T.tmp" "$T"
fi
run_wl noise3_nb2_LHC_B "$SCL/noise3_nb2_LHC_B.log" "$NK" "$WL/run_nb2_distributions.wls" "$T" LHC B 1000000 "$NK"
run_wl noise3_nb3_FASER2-FPF "$SCL/noise3_nb3_FASER2-FPF.log" "$NK" "$WL/run_nb3_sensitivity.wls" "$T" FASER2-FPF "3" "$NK"
SD="Sensitivity domains/ALP-fermion/FASER2-FPF"
V="$SCW/validation/noise_faser2.txt"; : > "$V"
M=0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.7,1.0,1.5,2.0,2.5,3.0,3.5
HAVE_CAMPAIGN=0
ref_or_skip "$V" "FASER2-FPF noise and production runs vs the campaign" "$CAMPAIGN_TREE" && HAVE_CAMPAIGN=1
for nev in "3." "2.3"; do
  F="Sensitivity_ALP-fermion_at_FASER2-FPF_Nev=${nev}_Npot=2.16e17.json"
  for pair in "$T/$SD/$F|$P/$SD/$F|noise run vs production run (both Linux)" \
              "$T/$SD/$F|$CAMPAIGN_TREE/$SD/$F|noise run vs campaign" \
              "$P/$SD/$F|$CAMPAIGN_TREE/$SD/$F|production run vs campaign"; do
    IFS='|' read -r a b l <<< "$pair"
    [ "$HAVE_CAMPAIGN" = 1 ] || [ "${b#"$CAMPAIGN_TREE"}" = "$b" ] || continue
    python3 -I "$PY/compare_curves.py" "$a" "$b" --variants central --masses "$M" --tol 0.01 --label "N>=$nev $l" \
      --out "$SCW/validation/noise_tmp.tsv" | grep -v '^  central within' | tee -a "$V" || true
    grep '^AP-Central|Central' "$SCW/validation/noise_tmp.tsv" | grep -v range | \
      awk -F'\t' '{l=($4=="-")?"   -  ":sprintf("%+.2f",($4-1)*100); u=($6=="-")?"   -  ":sprintf("%+.2f",($6-1)*100); printf "    m=%-5s lower %s%%  upper %s%%\n", $2, l, u}' | tee -a "$V"
  done
done
rm -f "$SCW/validation/noise_tmp.tsv"
