#!/usr/bin/env bash
# One BC4 or BC5 scan job: one mass, one worker, the published settings (scalar_eventcalc_20261019/lxplus/job_sc.sh):
#   bc4: python -m grendel.models.scalar.scan --mass M --workers 1 --thresholds 3 10 --vectors-dir V \
#          --templates-dir T --geometry-dir G --out O
#   bc5: python -m grendel.models.scalar.scan_bc5 --workers 1 --decay-samples 100 --thresholds 3 10 --vectors-dir V \
#          --templates-dir T --geometry-dir G --out O --mass M
#
#   scan_one.sh <bc4|bc5> <mass> [RUN]        RUN (default: main) names the output tree work/scalar/<bench>/scans/<RUN>
#
# The vectors are work/scalar/<bench>/llp_4vectors; the templates are work/scalar/templates_exhad (ours, Linux), or
# for RUN=control_reftpl the published exHad bundles (macOS arm64, v2_baseline_20260929/templates/scalar of the
# reference data, $GRENDEL_REF; without it a control_reftpl job is skipped). Each mass writes its own folder
# <RUN>/out/<mass> (sensitivity.csv, run_metadata.json, scan.log, time.txt); a finished mass is skipped, an
# interrupted one restarts (--resume). The ray-cast cache is shared by the runs of a benchmark (one file per mass).
# Holds one slot of the shared pool; scan_job.sh starts the jobs one by one when 3 GB of memory are available.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
B=$1; M=$2; RUN=${3:-main}
T="$W/templates_exhad"
V="$W/$B/llp_4vectors"
G="$W/$B/geometry_cache"
O="$W/$B/scans/$RUN/out/$M"
[ -s "$O/sensitivity.csv" ] && [ -s "$O/time.txt" ] && { echo "$B $M ($RUN): done"; exit 0; }
if [ "$RUN" = control_reftpl ]; then   # the published bundles are reference data ($GRENDEL_REF)
  have_ref "$B $M control_reftpl scan (the published exHad bundles)" "$REF_TEMPLATES" || exit 0
  T="$REF_TEMPLATES"
fi
mkdir -p "$O" "$G"
case $B in
  bc4) CMD=(python -m grendel.models.scalar.scan --mass "$M" --workers 1 --thresholds 3 10 --vectors-dir "$V"
            --templates-dir "$T" --geometry-dir "$G" --out "$O" --resume) ;;
  bc5) CMD=(python -m grendel.models.scalar.scan_bc5 --workers 1 --decay-samples 100 --thresholds 3 10 --vectors-dir "$V"
            --templates-dir "$T" --geometry-dir "$G" --out "$O" --mass "$M" --resume) ;;
  *) echo "usage: scan_one.sh <bc4|bc5> <mass> [RUN]"; exit 2 ;;
esac
# one job per mass at a time, even if two dispatchers (07 and 08) ask for the same mass
exec 9>"$O/.lock"
flock 9
[ -s "$O/sensitivity.csv" ] && [ -s "$O/time.txt" ] && { echo "$B $M ($RUN): done"; exit 0; }
$SLOT -n 1 --label "scalar:scan_${B}_${M}_${RUN}" -- bash "$S/scan_job.sh" "$O" "${CMD[@]}" && rc=0 || rc=$?
[ $rc -eq 0 ] && [ -s "$O/sensitivity.csv" ] && echo "$B $M ($RUN): $(cat "$O/time.txt")" || echo "$B $M ($RUN): FAILED rc=$rc"
exit $rc
