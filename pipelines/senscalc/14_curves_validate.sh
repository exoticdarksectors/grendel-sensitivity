#!/usr/bin/env bash
# 14_curves_validate.sh -- curve files in the format curves/ reads, and validation against the campaign. No kernel needed.
#   source env/activate.sh && bash pipelines/senscalc/14_curves_validate.sh
# Inputs : work/senscalc/run/SensCalc@eb2c456/Sensitivity domains/ALP-fermion/{FASER2-FPF,FASER2-FPF-Etot100,SHiP-ECN3}/*.json
# Outputs: work/senscalc/out/faser2_fpf_senscalc_nev3.dat, _nev2p3.dat   (campaign make_curve_files.py, same header)
#          work/senscalc/out/projections/faser2_senscalc_fpf.dat     -> would replace curves/data/bc10/projections/faser2_senscalc_fpf.dat
#          work/senscalc/out/projections/ship_senscalc.dat           -> regenerated from M. Ovchynnikov's JSON (must be byte-identical);
#                                                                       that JSON is not public (reference tree only): without
#                                                                       it, a copy of the tracked file
#          work/senscalc/out/projections/ship_senscalc_rerun.dat     -> our SHiP rerun (validation only)
#          work/senscalc/validation/*.tsv, summary.txt
# The comparisons with the campaign and with M. Ovchynnikov's JSON need the optional reference tree ($GRENDEL_REF);
# without it they are skipped (SKIP lines in summary.txt). Those with the tracked files always run.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
T="${1:-$SCW/run/SensCalc@eb2c456}"
OUT="${SENSCALC_OUT:-$SCW/out}"; VAL="${SENSCALC_VAL:-$SCW/validation}"; mkdir -p "$OUT/projections" "$VAL"
SD="Sensitivity domains/ALP-fermion"
TRACKED="$REPRO/src/grendel-sensitivity/curves/data/bc10/projections"
S="$VAL/summary.txt"; : > "$S"
say() { echo "$@" | tee -a "$S"; }

# ---- curve files ----
python3 -I "$PY/make_curve_files.py" "$T" "$OUT" | tee -a "$S"
cp "$OUT/faser2_fpf_senscalc_nev3.dat" "$OUT/projections/faser2_senscalc_fpf.dat"
if ref_or_skip "$S" "SHiP curve ship_senscalc.dat from M. Ovchynnikov's JSON (not public)" "$MAKSYM_JSON"; then
  MJ=$MAKSYM_JSON
else
  MJ=-
fi
( cd "$REPRO/src/grendel-sensitivity" && python -I "$PY/ship_curves.py" "$OUT" "$MJ" \
    "$T/$SD/SHiP-ECN3/Sensitivity_ALP-fermion_at_SHiP-ECN3_Nev=2.3_Npot=6.e20.json" ) | tee -a "$S"
if [ "$MJ" = - ]; then
  cp "$TRACKED/ship_senscalc.dat" "$OUT/projections/ship_senscalc.dat"
  say "ship_senscalc.dat: copy of the tracked curves/data/bc10/projections/ship_senscalc.dat (not regenerated)"
fi

say "== byte checks"
for f in ship_senscalc.dat faser2_senscalc_fpf.dat; do
  [ "$MJ" = - ] && [ "$f" = ship_senscalc.dat ] && continue
  a=$(sha256sum < "$OUT/projections/$f" | cut -c1-16); b=$(sha256sum < "$TRACKED/$f" | cut -c1-16)
  say "$f: ours $a tracked $b $([ "$a" = "$b" ] && echo IDENTICAL || echo DIFFERENT)"
done

say "== FASER2-FPF and SHiP-ECN3 regions vs the campaign (all 9 variants; MC noise of the campaign's own rerun: <= 0.5 %)"
cmp() { python3 -I "$PY/compare_curves.py" "$@" | tee -a "$S" || true; }
if ref_or_skip "$S" "FASER2-FPF curves vs the campaign" "$CAMPAIGN_TREE"; then
  for e in FASER2-FPF FASER2-FPF-Etot100; do
    for nev in "3." "2.3"; do
      f="Sensitivity_ALP-fermion_at_${e}_Nev=${nev}_Npot=2.16e17.json"
      cmp "$T/$SD/$e/$f" "$CAMPAIGN_TREE/$SD/$e/$f" --label "$e N>=$nev: ours vs campaign" --tol 0.01 --out "$VAL/${e}_nev${nev%.}_vs_campaign.tsv"
    done
  done
fi
cmp "$OUT/projections/faser2_senscalc_fpf.dat" "$TRACKED/faser2_senscalc_fpf.dat" --label "faser2_senscalc_fpf.dat: ours vs tracked" \
  --tol 0.01 --out "$VAL/faser2_senscalc_fpf_dat_vs_tracked.tsv" \
  --masses 0.03,0.05,0.1,0.2,0.3,0.5,0.7,1.0,1.2,1.5,1.8,2.0,2.5,3.0,3.3,3.5,3.6
f="Sensitivity_ALP-fermion_at_SHiP-ECN3_Nev=2.3_Npot=6.e20.json"
if ref_or_skip "$S" "SHiP-ECN3 rerun vs M. Ovchynnikov's JSON" "$MAKSYM_JSON"; then
  cmp "$T/$SD/SHiP-ECN3/$f" "$MAKSYM_JSON" --label "SHiP-ECN3 rerun vs M. Ovchynnikov 2026-10-16" --tol 0.005 --out "$VAL/ship_rerun_vs_maksym.tsv"
fi
if ref_or_skip "$S" "SHiP-ECN3 rerun vs the campaign rerun" "$CAMPAIGN_TREE"; then
  cmp "$T/$SD/SHiP-ECN3/$f" "$CAMPAIGN_TREE/$SD/SHiP-ECN3/$f" --label "SHiP-ECN3 rerun vs campaign rerun" --tol 0.01 --out "$VAL/ship_rerun_vs_campaign.tsv"
fi
if ref_or_skip "$S" "(reference) campaign SHiP-ECN3 rerun vs M. Ovchynnikov's JSON" "$CAMPAIGN_TREE" "$MAKSYM_JSON"; then
  cmp "$CAMPAIGN_TREE/$SD/SHiP-ECN3/$f" "$MAKSYM_JSON" --label "(reference) campaign rerun vs M. Ovchynnikov" --tol 0.005 --out "$VAL/ship_campaign_vs_maksym.tsv"
fi

say "== deterministic channels (bremsstrahlung, fragmentation) vs the campaign, tabulated N_ev (py/smoke_check.py)"
if ref_or_skip "$S" "deterministic channels vs the campaign" "$CAMPAIGN_TREE"; then
  for e in FASER2-FPF FASER2-FPF-Etot100 SHiP-ECN3; do
    say "-- $e"; python3 -I "$PY/smoke_check.py" "$T" "$CAMPAIGN_TREE" "$e" | grep -v '^curve\|^SMOKE' | tee -a "$S" || true
  done
fi

say "== deterministic metadata of the acceptance tables (Br_vis, ctau(m, y=1), particle numbers) vs the campaign"
if ref_or_skip "$S" "acceptance metadata vs the campaign" "$CAMPAIGN/logs/inspect_acceptances.txt"; then
  for e in SHiP-ECN3 FASER2-FPF; do
    [ -s "$SCL/inspect_acceptance_$e.txt" ] || continue
    python3 -I "$PY/compare_inspect.py" "$SCL/inspect_acceptance_$e.txt" "$CAMPAIGN/logs/inspect_acceptances.txt" "$e" | tee -a "$S"
  done
fi

say "== channel breakdown at the lower edge (FASER2-FPF, N>=3)"
if python3 -I "$PY/channel_breakdown.py" "$T/Tabulated Nevents/ALP-fermion/FASER2-FPF" 3 0.1 0.3 0.5 1 1.5 2 3 3.5 3.7 \
     > "$VAL/channel_breakdown_FASER2-FPF.txt" \
   && ref_or_skip "$S" "channel breakdown vs the campaign" "$CAMPAIGN/logs/channel_breakdown_FASER2-FPF.txt"; then
  diff <(grep '^m=' "$VAL/channel_breakdown_FASER2-FPF.txt" | cut -c1-60) \
    <(grep '^m=' "$CAMPAIGN/logs/channel_breakdown_FASER2-FPF.txt" | cut -c1-60) | tee -a "$S" || true
fi
echo "summary: $S"
