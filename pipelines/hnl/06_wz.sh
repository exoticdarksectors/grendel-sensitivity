#!/bin/bash
# HNL step 06: MadGraph W/Z -> l N per (flavour, mass) on copies of the step-03 process directories (lib/run_wz_point.py):
#   coarse  312 points, 100 000 events, seed offset r=21: the June v1 W/Z (central_newgrids_20260623/wz_fix.log,
#           MG5 3.6.6, --nb_core=1, a freshly generated directory per mass, i.e. the template randinit r=21)
#   fine    411 points, 1 000 000 events, r=27 (Ue, Utau), r=33 (Umu), r=34 (Umu 3.845): the lxplus runs of
#           hnl_fine_condor_20261014 (results/raw/*/generate_events_run_*.log), gzipped (1M rows each)
# K = 1.3 on every row (run_wz_production). One slot per job, --nb_core=1: MadGraph gives the same events for any
# nb_core (checked at Umu 3 GeV, 1M: nb_core 1 and 4 identical). Resume-safe.   06_wz.sh [coarse|fine|both] [jobs]
set -euo pipefail
source "$(dirname "$0")/../../env/activate.sh"
G=${1:-both}; N=${2:-2}
grids="coarse fine"; [ $G = both ] || grids=$G
python $REPRO/pipelines/hnl/lib/run_points.py wz --grid $grids --jobs $N --slots 1
