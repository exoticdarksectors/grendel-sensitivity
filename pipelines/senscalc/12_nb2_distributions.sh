#!/usr/bin/env bash
# 12_nb2_distributions.sh -- SensCalc notebook 2 (ALP angle-energy distributions), headless, 10^6 events per channel.
#   source env/activate.sh && bash pipelines/senscalc/12_nb2_distributions.sh
# Reproduces the campaign's nb2 runs: facility LHC (FASER2-FPF) and SPS (SHiP-ECN3), channels
# B, Eta, EtaPr, KS, Omega, RhoCh, Bremsstrahlung(-AP); Fragmentation and Drell-Yan are pregenerated in the repository.
# Output: work/senscalc/run/SensCalc@eb2c456/spectra/New physics particles spectra/ALP-fermion/DoubleDistr_ALP-fermion_<mode>_<facility>.m
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
require_kernel
T="$SCW/run/SensCalc@eb2c456"
make_tree "$T"
MODES=B,Eta,EtaPr,KS,Omega,RhoCh,Bremsstrahlung
run_wl nb2_LHC "$SCL/nb2_LHC.log" "$NK" "$WL/run_nb2_distributions.wls" "$T" LHC "$MODES" 1000000 "$NK"
run_wl nb2_SPS "$SCL/nb2_SPS.log" "$NK" "$WL/run_nb2_distributions.wls" "$T" SPS "$MODES" 1000000 "$NK"
