#!/bin/bash
# HNL step 01: rebuild the charged-kaon SoftQCD spectrum (Pythia 8.317, 60 000 events, seed 42: the recipe of the
# committed grendel/models/hnl/production/data/kaon_softqcd_spectrum.npz) into work/hnl/kaon_spectrum/ and compare
# it with the committed file. The production (steps 03-05) reads the committed file, as the published production did.
# Resume-safe: skips when the rebuilt file exists.
set -euo pipefail
source "$(dirname "$0")/../../env/activate.sh"
O=$REPRO/work/hnl/kaon_spectrum; mkdir -p $O
L=$REPRO/logs/hnl/01_kaon_spectrum.log
[ -s $O/kaon_softqcd_spectrum.npz ] && { echo "exists: $O/kaon_softqcd_spectrum.npz"; exit 0; }
export CXX=${CXX:-c++}
mkdir -p $O/bin && ln -sf "$CXX" $O/bin/c++          # make_kaon_spectrum calls "c++": the conda g++ that built Pythia
cd $O
PATH=$O/bin:$PATH python $REPRO/pipelines/common/slot_run.py -n 1 --label hnl:01_kaon_spectrum -- \
  python -m grendel.models.hnl.production.make_kaon_spectrum -o $O/kaon_softqcd_spectrum.npz > $L 2>&1
python - $O/kaon_softqcd_spectrum.npz $REPRO/src/grendel-sensitivity/grendel/models/hnl/production/data/kaon_softqcd_spectrum.npz <<'PY' | tee -a $L
import sys, numpy as np
new, ref = np.load(sys.argv[1]), np.load(sys.argv[2])
for k in ref.files:
    a, b = new[k], ref[k]
    if k == "hist":
        d = np.count_nonzero(a != b)
        print(f"hist: {d}/{a.size} bins differ; sum {a.sum():.0f} vs {b.sum():.0f}")
    else:
        print(f"{k}: {'identical' if np.array_equal(a, b) else f'rebuilt {a!r} committed {b!r}'}")
PY
