#!/usr/bin/env bash
# 11_nb1_acceptances.sh -- SensCalc notebook 1 (acceptances, ALP-fermion, "Default" decays), headless, full statistics.
#   source env/activate.sh && bash pipelines/senscalc/11_nb1_acceptances.sh
# Reproduces campaigns/bc10/faser2_senscalc_20261025 (README "Drivers"):
#   SHiP-ECN3           cuts "ship"    (docs/recommended-selections.md: E > 1 GeV every category, IP < 2.5 m)  [7 min on M-series]
#   FASER2-FPF          cuts "none"    (">=2" detectable products, no kinematic cut)                         [3 min]
#   FASER2-FPF-Etot100  cuts "etot100" (E of the detected pair > 100 GeV; selection variant)                [5 min]
# Output: work/senscalc/run/SensCalc@eb2c456/Acceptances/<name>/Acceptance_<name>_for_ALP-fermion_decayDescription_Default.m
# Resume-safe (markers in work/senscalc/markers).
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
require_kernel
T="$SCW/run/SensCalc@eb2c456"
make_tree "$T"
run_wl nb1_SHiP-ECN3_ship "$SCL/nb1_SHiP-ECN3_ship.log" "$NK" "$WL/run_nb1_acceptance.wls" "$T" SHiP-ECN3 SHiP-ECN3 ship "$NK"
run_wl nb1_FASER2-FPF_none "$SCL/nb1_FASER2-FPF_none.log" "$NK" "$WL/run_nb1_acceptance.wls" "$T" FASER2-FPF FASER2-FPF none "$NK"
run_wl nb1_FASER2-FPF_etot100 "$SCL/nb1_FASER2-FPF_etot100.log" "$NK" "$WL/run_nb1_acceptance.wls" "$T" FASER2-FPF FASER2-FPF-Etot100 etot100 "$NK"
# metadata of the three acceptance tables (facility, ECAL, mass grid, Br_vis, ctau(m, y=1))
for e in SHiP-ECN3 FASER2-FPF FASER2-FPF-Etot100; do
  f="$T/Acceptances/$e/Acceptance_${e}_for_ALP-fermion_decayDescription_Default.m"
  [ -s "$SCL/inspect_acceptance_$e.txt" ] || with_slots 1 inspect wolframscript -file "$WL/inspect_acceptance.wls" "$f" > "$SCL/inspect_acceptance_$e.txt" 2>&1 < /dev/null
done
