#!/usr/bin/env bash
# Step 08 -- the full-statistics production: every BC4 (86) and BC5 (112) mass, one job per mass (scan_one.sh,
# RUN=main), with our four-vectors (01, 05) and Linux exHad templates (06).
#
# Published recipe: scalar_eventcalc_20261019 (HTCondor clusters 16501879/16501880, lxplus/job_sc.sh): code c621279
# (= the grendel/ of 8802f50 for these scans), one worker, thresholds 3 10, BC5 --decay-samples 100, per-mass seeds.
# lxplus took 31.7 (BC4) + 39.7 (BC5) core-hours; see the README for the measured times here.
# Jobs are dispatched SCALAR_JOBS (default 6) at a time, BC5 first (the longer ones), each waiting for one slot of the
# shared pool and for 3 GB of free memory. Finished masses are skipped, so rerunning this script resumes.
# Output: work/scalar/<bench>/scans/main/out/<mass>/; log logs/scalar/08_scan_all.log.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
step_lock
J=${SCALAR_JOBS:-6}
{
  for m in $BC5_MASSES; do echo "bc5 $m main"; done
  for m in $BC4_MASSES; do echo "bc4 $m main"; done
} | xargs -P "$J" -L 1 bash "$S/scan_one.sh" 2>&1 | tee -a "$L/08_scan_all.log"
n4=$(find "$W/bc4/scans/main/out" -name time.txt | wc -l); n5=$(find "$W/bc5/scans/main/out" -name time.txt | wc -l)
echo "$(date '+%F %T') 08_scan_all: bc4 $n4/86, bc5 $n5/112" | tee -a "$L/08_scan_all.log"
[ "$n4" -eq 86 ] && [ "$n5" -eq 112 ] || { echo "08_scan_all INCOMPLETE (rerun to resume)"; exit 1; }
echo "08_scan_all OK"
