#!/usr/bin/env bash
# Step 03: the NNPDF4.0 raw quark grids (fonllgridlha) of the central grids, 100 rapidities x 80 pT nodes
# (0.05-200 GeV) per quark, split into FONLL_NCHUNKS (default 4) contiguous rapidity chunks of 25 rows, each its
# own process and slot; FONLL_PARALLEL (default 4) chunks at a time. ~1.5-2.5 h per chunk here, ~6-9 single-core
# hours per quark.
#
#   source env/activate.sh && bash pipelines/fonll/03_raw_quark_grids.sh [bottom] [charm]
#
# Card: NNPDF40_nlo_as_01180 (331700), pp 7000 + 7000 GeV, m_b = 4.75 / m_c = 1.50 GeV, Lambda_5 default,
# ffact = fren = 1, the 100 rapidities of generate.Y_VALUES, pT_max 200 GeV with 80 automatic nodes, FONLL
# (icalctype 1): gen.grid_input, as `generate --pdf nlo --quark <q> --grid-workers 4` (raw_chunk.py says why the
# chunks are separate processes). Output: work/fonll/grids/run/nlo_<q>_central/nl_<q>_.out (chunks in
# grid_chunks_n4/), compared with fonll-local/run/nlo_<q>/nl_<q>_.out (the raw grids of the tracked grids).
set -euo pipefail
: "${REPRO:?source env/activate.sh first}"
QUARKS=("$@"); [ ${#QUARKS[@]} -gt 0 ] || QUARKS=(bottom charm)
NCH="${FONLL_NCHUNKS:-4}"
LOG="$REPRO/logs/fonll/03_raw_quark_grids.log"
mkdir -p "$(dirname "$LOG")"
exec > >(tee -a "$LOG") 2>&1
echo "=== $(date -Is) 03_raw_quark_grids: ${QUARKS[*]}, $NCH chunks, ${FONLL_PARALLEL:-4} at a time"
T="$REPRO/pipelines/fonll/task.sh"
for q in "${QUARKS[@]}"; do for ((i = 0; i < NCH; i++)); do echo "raw:$q:$i"; done; done \
  | xargs -P "${FONLL_PARALLEL:-4}" -n 1 bash "$T"
for q in "${QUARKS[@]}"; do bash "$T" "merge:$q"; done
echo "=== $(date -Is) 03_raw_quark_grids done"
