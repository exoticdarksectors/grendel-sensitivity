#!/usr/bin/env bash
# 10_smoke.sh -- smoke run of the headless SensCalc chain at tiny statistics (FASER2-FPF only).
#   source env/activate.sh && bash pipelines/senscalc/10_smoke.sh
# Tree: work/senscalc/smoke/SensCalc@eb2c456 (separate from the production tree).
#   nb1  acceptance FASER2-FPF, ">=2", no cuts (notebook-default statistics; the acceptance has no size knob in the driver)
#   nb2  LHC distributions with ${SMOKE_NSIM:-20000} events per channel (production: 10^6)
#   nb3  FASER2-FPF sensitivity, N_ev 2.3 and 3
# Then checks that the chain ran and landed in the right place (py/smoke_check.py; needs the campaign of the
# optional reference tree, skipped without it).
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
require_kernel
T="$SCW/smoke/SensCalc@eb2c456"
NSIM="${SMOKE_NSIM:-20000}"
make_tree "$T"
run_wl smoke_nb1_FASER2-FPF "$SCL/smoke_nb1_FASER2-FPF.log" "$NK" "$WL/run_nb1_acceptance.wls" "$T" FASER2-FPF FASER2-FPF none "$NK"
run_wl smoke_nb2_LHC "$SCL/smoke_nb2_LHC.log" "$NK" "$WL/run_nb2_distributions.wls" "$T" LHC B,Eta,EtaPr,KS,Omega,RhoCh,Bremsstrahlung "$NSIM" "$NK"
run_wl smoke_nb3_FASER2-FPF "$SCL/smoke_nb3_FASER2-FPF.log" "$NK" "$WL/run_nb3_sensitivity.wls" "$T" FASER2-FPF "3" "$NK"
# criterion (py/smoke_check.py): deterministic channels exact (bremsstrahlung 0.5 %, fragmentation 2 %), curve placement
# sized for 2e4 events (edges within 50 %, reach within 5 %); a percent-level curve tolerance is meaningless at this
# statistic because FASER2 sees only a few sampled B-decay ALPs per mass (see the docstring there)
if have_ref "SensCalc smoke run vs the campaign (deterministic channels, curve placement)" "$CAMPAIGN_TREE"; then
  python3 -I "$PY/smoke_check.py" "$T" "$CAMPAIGN_TREE" FASER2-FPF | tee "$SCL/10_smoke_compare.log"
else
  echo "smoke run finished (all three notebook runs done); no criterion without the campaign"
fi
