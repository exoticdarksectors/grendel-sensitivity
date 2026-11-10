#!/bin/bash
# HNL step 08: one scan job per (flavour, mass) (lib/scan_point.py; the lxplus jobs of v2_3_20261016).
# A point waits (holding no slot) until its inputs exist; fine points also until lib/memgate.py lets them start.
#   08_scans.sh [coarse|fine|both] [concurrent jobs] [more lib/run_points.py options, e.g. --give-up-when FILE...]
set -euo pipefail
source "$(dirname "$0")/../../env/activate.sh"
G=${1:-both}; N=${2:-2}; shift $(( $# < 2 ? $# : 2 ))
grids="coarse fine"; [ $G = both ] || grids=$G
python $REPRO/pipelines/hnl/lib/run_points.py scan --grid $grids --jobs $N "$@"
