#!/usr/bin/env bash
# Step 06 -- Fig. 7 (BC10): render the paper panel from our table and from the published table, and compare both
# PNGs pixel by pixel with the paper's figure ($GRENDEL_PAPER_FIGURES/bc10_grendel_paper.png, in this repository).
#   python -m curves.plot paper --grendel-dir <tables> --benchmark bc10 --out <dir>
# The comparison curves (SHiP, FASER2, bounds) are the tracked files of src/grendel-sensitivity/curves/data.
#   06_plot.sh [VARIANT]     (default main: work/alp/tables_main)
# Outputs: work/alp/figures_<VARIANT>/, work/alp/figures_published/, work/alp/validation/figure_<VARIANT>.txt.
# Seconds; no slot needed.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
VARIANT="${1:-main}"
cd "$W"
python -m curves.plot paper --grendel-dir "$W/tables_$VARIANT" --benchmark bc10 --out "$W/figures_$VARIANT" \
  > "$L/06_plot_$VARIANT.log" 2>&1
python -m curves.plot paper --grendel-dir "$(dirname "$PUB")" --benchmark bc10 --out "$W/figures_published" \
  > "$L/06_plot_published.log" 2>&1
python - "$PAPER_PNG" "$W/figures_published/bc10_grendel_paper.png" "$W/figures_$VARIANT/bc10_grendel_paper.png" \
  > "$W/validation/figure_$VARIANT.txt" <<'EOF'
import hashlib, sys
import numpy as np
from PIL import Image
import os
paper, *others = sys.argv[1:]
def rel(p):  # paths relative to the checkout, so the report names no machine-specific folder
    return os.path.relpath(p, os.environ["REPRO"])
ref = np.asarray(Image.open(paper).convert("RGB"), dtype=int)
print(f"paper PNG {rel(paper)}: {ref.shape[1]}x{ref.shape[0]}, sha256 {hashlib.sha256(open(paper, 'rb').read()).hexdigest()[:16]}")
for p in others:
    img = np.asarray(Image.open(p).convert("RGB"), dtype=int)
    if img.shape != ref.shape:
        print(f"{rel(p)}: size {img.shape} differs from the paper PNG {ref.shape}")
        continue
    d = np.abs(img - ref).max(axis=2)
    print(f"{rel(p)}: sha256 {hashlib.sha256(open(p, 'rb').read()).hexdigest()[:16]}; pixels differing from the paper PNG "
          f"{int((d > 0).sum())} of {d.size} ({(d > 0).mean():.2e}), > 32 levels {int((d > 32).sum())}, max {d.max()}")
EOF
cat "$W/validation/figure_$VARIANT.txt"
