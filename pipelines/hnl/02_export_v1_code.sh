#!/bin/bash
# HNL step 02: export the code of the June 2026 v1 HNL production, whose B-meson, Lambda_b, prompt-tau and induced-tau
# channel files the published coarse-grid rows still carry ("the other channels of the v2 production",
# curves/data/grendel_v2_3/README.md; asbuilt_20261014/hnl/llp_4vectors/*/{Bmeson,Bbaryon,tau,induced_tau} link to
# signal-models/hnl/llpatcolliders_HNL_PR1/hnl/tmp/runs/central_newgrids_20260623, produced 2026-06-23 ~12:30).
# The code at that time is llpatcolliders 06e5a18 (2026-06-23 10:38, branch hnl-production-mesons-ew, the last HNL
# commit before the run): git archive of its hnl/ into src/llpatcolliders-hnl-06e5a18. It vendors HNLCalc
# leoredi/HNLCalc e292cef9 (G_F = 1.166378e-5; the pin b1d36c6 has 1.1663788e-5 and the FLAG 2024 f_D, f_Ds) and
# has SIGMA_BC_PB = 0.9e6 (the B_c -> tau nu source of the induced-tau pool; 1.64e6 since v2).
#
# The export is checked in under src/, so this step only checks it (COMMIT file). Without it, it re-creates the
# export from the public repository https://github.com/leoredi/llpatcolliders (branch hnl-production-mesons-ew
# contains 06e5a18). Resume-safe.
#     02_export_v1_code.sh [destination]      (default src/llpatcolliders-hnl-06e5a18)
set -euo pipefail
source "$(dirname "$0")/../../env/activate.sh"
URL=https://github.com/leoredi/llpatcolliders
C=06e5a184b348ec70923e8a5ecf75579a05e3a057
D=${1:-$REPRO/src/llpatcolliders-hnl-06e5a18}
if [ -f "$D/hnl/run_all.py" ] && [ "$(head -1 "$D/COMMIT" 2>/dev/null)" = $C ]; then
  echo "exists: $D ($(sed -n 2p "$D/COMMIT"))"; exit 0
fi
if [ -e "$D" ] && [ -n "$(ls -A "$D")" ]; then
  echo "$D exists but is not the export of $C (no hnl/run_all.py or another COMMIT): move it away first" >&2; exit 1
fi
G=$(mktemp -d); trap 'rm -rf "$G"' EXIT
git init -q "$G"
git -C "$G" fetch -q --depth 1 "$URL" $C 2>/dev/null || git -C "$G" fetch -q "$URL" hnl-production-mesons-ew
git -C "$G" cat-file -e "$C^{commit}"
mkdir -p "$D"
git -C "$G" archive --format=tar $C hnl | tar -x -C "$D"
{ git -C "$G" rev-parse $C; git -C "$G" log -1 --format='%H %ad %s' --date=iso $C; } > "$D/COMMIT"
echo "exported $(head -1 "$D/COMMIT") -> $D"
