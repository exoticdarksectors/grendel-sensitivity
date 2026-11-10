# Sourced by every pipelines/fonll step (after env/activate.sh): where the FONLL area lives.
: "${REPRO:?source env/activate.sh first}"
export FONLL_WORK="$REPRO/work/fonll"
# the patched FONLL tree built by 01_build_fonll.sh (third_party/fonll stays pristine)
export GRENDEL_FONLL_DIR="${GRENDEL_FONLL_DIR:-$FONLL_WORK/fonll-1.3.3}"
# central grids, run scratch and per-command logs of generate.py
export GRENDEL_FONLL_GRID_DIR="$FONLL_WORK/grids/output"
export FONLL_RUN_DIR="$FONLL_WORK/grids/run"
export FONLL_LOG_DIR="$FONLL_WORK/grids/logs"
# NNPDF40_nlo_as_01180 from the shared work/lhapdf, cteq66 from 00_fetch_pdfs.sh. generate.py copies
# GRENDEL_LHAPDF_DATA into LHAPDF_DATA_PATH, so both carry the two directories.
export LHAPDF_DATA_PATH="$REPRO/work/lhapdf:$FONLL_WORK/lhapdf"
export GRENDEL_LHAPDF_DATA="$LHAPDF_DATA_PATH"
# reference material (read-only): the tracked grids (always compared) and the original macOS FONLL workspace, in the
# optional reference tree (pipelines/common/ref.sh); the comparisons with the latter are skipped without it
source "$REPRO/pipelines/common/ref.sh"
export FONLL_TRACKED_DIR="$REPRO/src/grendel-sensitivity/grendel/production/data/fonll/central"
export FONLL_ORIG_WS="${FONLL_ORIG_WS:-$REF/grendel/shared/NNPDF40/fonll-local}"
mkdir -p "$REPRO/logs/fonll"
