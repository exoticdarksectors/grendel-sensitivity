#!/usr/bin/env bash
# Step 01: patch and build FONLL 1.3.3 (fonllgridlha + fragmfonll) with the conda gfortran and LHAPDF 6.5.6.
#
#   source env/activate.sh && bash pipelines/fonll/01_build_fonll.sh
#
# The pristine checkout third_party/fonll (alisw/fonll c7086e4) is copied to work/fonll/fonll-1.3.3 so that
# third_party stays pristine; the copy is patched and built by the shipped installer
#     python -m grendel.production.fonll_grids.install --fonll work/fonll/fonll-1.3.3 --build
# which checks the SHA-256 of the four files before and after patching against patches/PINS.json.
#
# Build flags. On Linux the shipped installer passes F77="gfortran -Wall -fno-automatic -fallow-argument-mismatch"
# and SYSOBJ= to make, but the Linux block of misc1/Makefile leaves OPT unset (-O0) and sets DEBUG=-ggdb, while
# the macOS build that produced the tracked grids compiled with OPT=-O (Darwin block). FONLL_OPT (default "-O")
# is therefore passed as OPT= on the make command line so the optimisation level matches the original build;
# FONLL_OPT="" reproduces the bare installer (-O0). The installer is called through a small wrapper that only
# appends OPT=$FONLL_OPT to its make command (install.MAKE); the code in src/ is not modified.
#
# Resume-safe: an existing build with the same flags is kept (make is incremental); a flag change rebuilds.
set -euo pipefail
: "${REPRO:?source env/activate.sh first}"
cd "$REPRO"

FONLL_SRC="$GRENDEL_THIRD_PARTY_DIR/fonll"
FONLL_DIR="${FONLL_BUILD_DIR:-$REPRO/work/fonll/fonll-1.3.3}"
FONLL_OPT="${FONLL_OPT--O}"
LOG="$REPRO/logs/fonll/01_build_fonll.log"
mkdir -p "$(dirname "$LOG")" "$REPRO/work/fonll"

exec > >(tee -a "$LOG") 2>&1
echo "=== $(date -Is) 01_build_fonll: OPT='$FONLL_OPT'"

if [ ! -d "$FONLL_DIR" ]; then
  echo "copy $FONLL_SRC -> $FONLL_DIR"
  mkdir -p "$(dirname "$FONLL_DIR")"
  rsync -a --exclude .git "$FONLL_SRC/" "$FONLL_DIR/"
fi

# a flag change invalidates the objects
STAMP="$FONLL_DIR/Linux/.grendel_build_flags"
WANT="F77=gfortran -Wall -fno-automatic -fallow-argument-mismatch | SYSOBJ= | OPT=$FONLL_OPT | $(gfortran --version | head -1) | LHAPDF $(lhapdf-config --version)"
if [ -f "$STAMP" ] && [ "$(cat "$STAMP")" != "$WANT" ]; then
  echo "build flags changed: cleaning $FONLL_DIR/Linux"
  rm -f "$FONLL_DIR"/Linux/*.o "$FONLL_DIR/Linux/fonllgridlha" "$FONLL_DIR/Linux/fragmfonll"
fi

python - "$FONLL_DIR" "$FONLL_OPT" <<'EOF'
import sys
from pathlib import Path
from grendel.production.fonll_grids import install
tree, opt = Path(sys.argv[1]), sys.argv[2]
install.MAKE += [f"OPT={opt}"]
print("make command:", " ".join(install.MAKE))
install.apply_patches(tree)
install.build(tree)
for rel, state in install.status(tree).items():
    print(f"  {state:9s} {rel}")
    assert state == "patched", rel
EOF

echo "$WANT" > "$STAMP"
sha256sum "$FONLL_DIR/Linux/fonllgridlha" "$FONLL_DIR/Linux/fragmfonll"
echo "=== $(date -Is) 01_build_fonll done"
