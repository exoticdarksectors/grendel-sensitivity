#!/usr/bin/env bash
# Step 03 -- BC5 Higgs production, stage 2: Pythia 8.317 parton shower (ISR + FSR, no MPI, no hadronisation, Higgs
# stable) of every LHE file, seed 1; records the final Higgs (pT, y, E, px, py, pz).
#
# Published recipe: asbuilt_20261014/logs/bc5_higgs_shower.log (code 461b7e5):
#   python -m grendel.models.scalar.higgs_shower          (all modes of the manifest, seed 1)
# higgs_shower.cc is compiled with the conda c++ against third_party/pythia8317 (no Python bindings needed; the
# original also used the C++ program, not pythia8mc).
# Two showers:
#   (a) of our LHE (02):                   work/scalar/bc5/madgraph/higgs_<mode>.csv
#   (b) of the stored 2026-09-25 LHE files: work/scalar/bc5/shower_of_published_lhe/higgs_<mode>.csv, a check of
#       the shower alone against the published asbuilt_20261014/bc5/madgraph/higgs_<mode>.csv. Needs the reference
#       data ($GRENDEL_REF: the stored LHE files and their shower); skipped without it. Nothing downstream needs (b)
#       except the check of 04.
# Resume-safe (a mode whose CSV is in the shower manifest is skipped). 1 slot, ~5 min per pass.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
step_lock
MG="$W/bc5/madgraph"
RS="$W/bc5/shower_of_published_lhe"
mkdir -p "$W/validation"
DIRS=("$MG")
if ref_ok "(b) shower of the stored 2026-09-25 LHE files vs the published shower" "$REF_SHOWER/higgs_lhe_manifest.json" "$REF_MG5"; then
  DIRS+=("$RS")
  mkdir -p "$RS"
  # (b) a manifest naming the stored LHE files (read-only) with their published LO cross sections
  python - "$REF_SHOWER/higgs_lhe_manifest.json" "$REF_MG5" "$RS/higgs_lhe_manifest.json" <<'EOF'
import json, sys
from pathlib import Path
src, ref_mg5, dst = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
m = json.loads(src.read_text())
for mode, e in m.items():
    e["lhe"] = str(ref_mg5 / f"higgs_{mode}" / "Events" / "run_200000" / "unweighted_events.lhe.gz")
    assert Path(e["lhe"]).exists(), e["lhe"]
dst.write_text(json.dumps(m, indent=2) + "\n")
EOF
fi
todo() {  # modes of $1 not yet showered
  python - "$1" <<'EOF'
import json, sys
from pathlib import Path
w = Path(sys.argv[1]); s = w / "higgs_shower_manifest.json"
done = json.loads(s.read_text()) if s.exists() else {}
print(" ".join(m for m in json.loads((w / "higgs_lhe_manifest.json").read_text())
               if not (m in done and Path(done[m]["csv"]).exists())))
EOF
}
for D in "${DIRS[@]}"; do
  modes=$(todo "$D")
  [ -z "$modes" ] && continue
  echo "$(date '+%F %T') 03_bc5_higgs_shower: $D modes $modes" | tee -a "$L/03_bc5_higgs_shower.log"
  # shellcheck disable=SC2086
  $SLOT -n 1 --label scalar:03_shower -- \
    python -m grendel.models.scalar.higgs_shower --work-dir "$D" --modes $modes --seed "$SHOWER_SEED" \
    >> "$L/03_bc5_higgs_shower.log" 2>&1
done
# the shower alone: published LHE showered here vs the published shower (macOS); only with the reference data
[ "${#DIRS[@]}" -eq 2 ] && python - "$RS" "$REF_SHOWER" <<'EOF' | tee "$W/validation/shower_of_published_lhe.txt"
import hashlib, sys
from pathlib import Path
import numpy as np
ours, ref = Path(sys.argv[1]), Path(sys.argv[2])
print("shower of the published LHE (Linux) vs the published shower (macOS), per mode:")
for mode in ("ggf", "vbf", "wh", "zh", "tth"):
    a, b = ours / f"higgs_{mode}.csv", ref / f"higgs_{mode}.csv"
    same = hashlib.sha256(a.read_bytes()).hexdigest() == hashlib.sha256(b.read_bytes()).hexdigest()
    A, B = np.loadtxt(a, delimiter=","), np.loadtxt(b, delimiter=",")
    if A.shape != B.shape:
        print(f"  {mode}: rows {len(A)} vs {len(B)}"); continue
    rows = np.any(A != B, axis=1)
    rel = np.max(np.abs(A - B) / np.maximum(np.abs(B), 1e-12), axis=1)
    print(f"  {mode}: sha256 equal {same}; {rows.sum()} of {len(A)} rows differ; rows differing by > 1e-6 rel: "
          f"{(rel > 1e-6).sum()}; <pT> {A[:, 0].mean():.4f} vs {B[:, 0].mean():.4f} GeV")
EOF
step_ok 03_bc5_higgs_shower
