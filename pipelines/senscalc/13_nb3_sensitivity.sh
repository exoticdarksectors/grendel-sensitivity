#!/usr/bin/env bash
# 13_nb3_sensitivity.sh -- SensCalc notebook 3 (ALP-fermion sensitivity): tabulated N_ev and sensitivity regions.
#   source env/activate.sh && bash pipelines/senscalc/13_nb3_sensitivity.sh
# Production "Revised (2501.04525)", all channels of the facility, default exposure (SHiP-ECN3 6e20 PoT; FASER2-FPF
# N_coll = 2.16e17), threshold 2.3 (notebook default) plus 3 for FASER2-FPF. All 9 Bremsstrahlung x Drell-Yan variants
# are written to each JSON.
# Output: work/senscalc/run/SensCalc@eb2c456/{Tabulated Nevents,Sensitivity domains}/ALP-fermion/<exp>/
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
require_kernel
T="$SCW/run/SensCalc@eb2c456"
for m in nb1_SHiP-ECN3_ship nb1_FASER2-FPF_none nb1_FASER2-FPF_etot100 nb2_LHC nb2_SPS; do
  [ -f "$MARK/$m" ] || { echo "ERROR: $m has not run (11_, 12_ first)"; exit 1; }
done
run_wl nb3_SHiP-ECN3 "$SCL/nb3_SHiP-ECN3.log" "$NK" "$WL/run_nb3_sensitivity.wls" "$T" SHiP-ECN3 "" "$NK"
run_wl nb3_FASER2-FPF "$SCL/nb3_FASER2-FPF.log" "$NK" "$WL/run_nb3_sensitivity.wls" "$T" FASER2-FPF "3" "$NK"
run_wl nb3_FASER2-FPF-Etot100 "$SCL/nb3_FASER2-FPF-Etot100.log" "$NK" "$WL/run_nb3_sensitivity.wls" "$T" FASER2-FPF-Etot100 "3" "$NK"
