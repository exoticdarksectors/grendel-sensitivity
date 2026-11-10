#!/usr/bin/env bash
# 22_check_quartic.sh -- BC5 quartic-coupling tables grendel/models/scalar/data/senscalc_quartic. No kernel needed.
#   source env/activate.sh && bash pipelines/senscalc/22_check_quartic.sh
# These are upstream SensCalc data files (phenomenology/Scalar/Production probabilities/), copied verbatim
# (PROVENANCE.md there). SensCalc only reads them (3. Scalar sensitivity.nb, codes/EventCalc/ProductionProbabilities.nb):
# no notebook in the repository generates them, so there is no Mathematica step to rerun. Reproduction = fetch + hash.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
TR="$REPRO/src/grendel-sensitivity/grendel/models/scalar/data/senscalc_quartic"
V="$SCW/validation/senscalc_quartic.txt"; : > "$V"
st=0
for f in BplustoXsSStotal.dat BstoSS.dat BplustoKplusSS.dat; do
  t=$(sha256sum < "$TR/$f" | cut -d' ' -f1)
  for tree in SensCalc SensCalc-v1.4; do
    u=$(sha256sum < "$REPRO/third_party/$tree/phenomenology/Scalar/Production probabilities/$f" | cut -d' ' -f1)
    r=$([ "$t" = "$u" ] && echo IDENTICAL || { st=1; echo DIFFERENT; })
    echo "$f  tracked ${t:0:16}  $tree $(git -C "$REPRO/third_party/$tree" rev-parse --short HEAD) ${u:0:16}  $r" | tee -a "$V"
  done
done
t=$(sha256sum < "$TR/LICENSE" | cut -d' ' -f1); u=$(sha256sum < "$REPRO/third_party/SensCalc-v1.4/LICENSE" | cut -d' ' -f1)
echo "LICENSE  $([ "$t" = "$u" ] && echo IDENTICAL || echo DIFFERENT) to SensCalc v1.4" | tee -a "$V"
exit $st
