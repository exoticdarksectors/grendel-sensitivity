#!/usr/bin/env bash
# Step 00: the LHAPDF sets the central FONLL grids need, besides NNPDF40_nlo_as_01180 (already in work/lhapdf,
# installed by env/bootstrap.sh):
#   cteq66 (LHAPDF id 10550): the charm D* -> D0 feeddown weight is fitted to public FONLL CTEQ6.6 D0 points
#   (generate.ensure_charm_calibration_grid), and the bottom public-point closure (step 04) uses it too.
#
#   source env/activate.sh && bash pipelines/fonll/00_fetch_pdfs.sh
#
# Downloaded from the LHAPDF upstream into work/fonll/lhapdf (kept apart from the shared work/lhapdf); the FONLL
# steps run with LHAPDF_DATA_PATH=work/lhapdf:work/fonll/lhapdf (pipelines/fonll/env.sh).
# Checks: the SHA-256 of the files the FONLL runs read (cteq66.info, cteq66_0000.dat) equal those of the copies
# in the original FONLL workspace (grendel/shared/NNPDF40/fonll-local/env/share/LHAPDF/cteq66 of the reference tree;
# the hashes are written below, so the check needs no reference data), and NNPDF40_nlo_as_01180 {.info, _0000.dat}
# equal the original ones as well. Resume-safe (skips what is there).
set -euo pipefail
: "${REPRO:?source env/activate.sh first}"
source "$REPRO/pipelines/fonll/env.sh"

DL="$FONLL_WORK/downloads"
DEST="$FONLL_WORK/lhapdf"
mkdir -p "$DL" "$DEST"
LOG="$REPRO/logs/fonll/00_fetch_pdfs.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== $(date -Is) 00_fetch_pdfs"

if [ ! -f "$DEST/cteq66/cteq66.info" ]; then
  if [ ! -f "$DL/cteq66.tar.gz" ]; then
    curl -fsSL -o "$DL/cteq66.tar.gz.part" https://lhapdfsets.web.cern.ch/lhapdfsets/current/cteq66.tar.gz
    mv "$DL/cteq66.tar.gz.part" "$DL/cteq66.tar.gz"
  fi
  tar -xzf "$DL/cteq66.tar.gz" -C "$DEST"
fi

# SHA-256 of the files the original (macOS) FONLL workspace read (fonll-local/env/share/LHAPDF)
L="$REPRO/work/lhapdf/NNPDF40_nlo_as_01180"
sha256sum -c - <<SUMS
d9c5eb32fddda6ea97c1ae65d0b0eac9eeff912d950b25cf2482731503f0c5d3  $DEST/cteq66/cteq66.info
6945d6e7537f6d6a49d2df39afdf1e65e11983e20903402addf1e1baf7d2671e  $DEST/cteq66/cteq66_0000.dat
9362949bda8c0ae1d6ba37542ee31ddbae010369d126ac45ab68ce8b306212a1  $L/NNPDF40_nlo_as_01180.info
b9d44f01c1ea60d229d94d370686ab3423a434a97adc5c26911a2f5ee7689ae8  $L/NNPDF40_nlo_as_01180_0000.dat
SUMS
echo "=== $(date -Is) 00_fetch_pdfs done"
