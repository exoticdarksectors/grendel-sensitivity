#!/usr/bin/env bash
# Gather the four rebuilt Sec. IV figures into figures/<area>/ and the rebuilt tables, SensCalc curves and validation
# reports into results/, and compare each figure with the paper's PNG (reference/paper_figures).
# Run by ../run_all.sh after every pipeline has finished; can be rerun on its own at any time (seconds, no slot).
#
#   figures/scalar/bc4_grendel_paper.{png,pdf}     Fig. 4   (pipelines/scalar, step 10)
#   figures/scalar/bc5_grendel_paper.{png,pdf}     Fig. 5   (pipelines/scalar, step 10)
#   figures/hnl/hnlimits_grendel_paper.{png,pdf}   Fig. 6   (pipelines/hnl, step 11)
#   figures/alp/bc10_grendel_paper.{png,pdf}       Fig. 7   rendered here: rebuilt table (pipelines/alp) and
#                                                           rebuilt SHiP/FASER2 SensCalc curves (pipelines/senscalc)
#
# Fig. 7 is also rendered from the rebuilt table with the tracked SensCalc curves (work/figures/bc10_tracked_curves/,
# the same figure as pipelines/alp step 06) and from the published table with the rebuilt curves
# (work/figures/bc10_published_table/), so that the effect of the table and of the SensCalc rerun can be told apart.
# The report is figures/comparison.txt; the per-figure diff images go to work/figures/diff/.
#
#   results/tables/{bc4,bc5,bc10,hnl}/sensitivity{,_nsig10}.csv   the rebuilt tables (layout of curves/data/grendel_v2_*)
#   results/senscalc/{ship_senscalc,faser2_senscalc_fpf}.dat      the rebuilt SensCalc curves of Fig. 7
#   results/validation/<pipeline>/                                the validation reports of work/<pipeline>/validation,
#                                                                 with this checkout's path written as ".", the
#                                                                 reference tree's as $GRENDEL_REF
set -uo pipefail
REPRO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [ -z "${GRENDEL_EXHAD_PYTHON:-}" ]; then
  # shellcheck disable=SC1091
  source "$REPRO/env/activate.sh"
fi
cd "$REPRO"
PAPER="${GRENDEL_PAPER_FIGURES#"$REPRO"/}"     # the figures of the paper (reference/paper_figures)
SRC="$REPRO/src/grendel-sensitivity"
W="$REPRO/work/figures"
REPORT="$REPRO/figures/comparison.txt"
CMP="$REPRO/pipelines/scalar/py/compare_png.py"
export MPLCONFIGDIR="$W/.mpl" PYTHONDONTWRITEBYTECODE=1
mkdir -p "$W/diff" figures/alp

rc=0
need() { if [ ! -s "$1" ]; then echo "MISSING $1 (from $2)" >&2; rc=1; return 1; fi; }

# --- Fig. 7: render from the rebuilt table, with the rebuilt and with the tracked SensCalc curves
TABLES="$REPRO/work/alp/tables_main"
SC="$REPRO/work/senscalc/out/projections"
if need "$TABLES/bc10/sensitivity.csv" "pipelines/alp step 05"; then
  ( cd "$SRC" && python -m curves.plot paper --grendel-dir "$TABLES" --benchmark bc10 \
      --out "$W/bc10_tracked_curves" ) > "$W/bc10_tracked_curves.log" 2>&1 || { echo "FAIL Fig. 7 (tracked curves)"; rc=1; }
  if [ -s "$SC/ship_senscalc.dat" ] && [ -s "$SC/faser2_senscalc_fpf.dat" ]; then
    ( cd "$SRC" && python "$REPRO/figures/render_bc10.py" "$TABLES" "$W/bc10_rebuilt_curves" \
        --ship "$SC/ship_senscalc.dat" --faser2 "$SC/faser2_senscalc_fpf.dat" ) > "$W/bc10_rebuilt_curves.log" 2>&1 \
      && cp "$W/bc10_rebuilt_curves/bc10_grendel_paper."{png,pdf} figures/alp/ \
      || { echo "FAIL Fig. 7 (rebuilt curves)"; rc=1; }
  else
    echo "NOTE pipelines/senscalc has not run (no Wolfram Engine?): Fig. 7 drawn with the tracked SensCalc curves" >&2
    cp "$W/bc10_tracked_curves/bc10_grendel_paper."{png,pdf} figures/alp/ || rc=1
  fi
fi
# the published table with the rebuilt SensCalc curves: the effect of the SensCalc rerun alone
if [ -s "$SC/ship_senscalc.dat" ] && [ -s "$SC/faser2_senscalc_fpf.dat" ]; then
  ( cd "$SRC" && python "$REPRO/figures/render_bc10.py" "$SRC/curves/data/grendel_v2_4" "$W/bc10_published_table" \
      --ship "$SC/ship_senscalc.dat" --faser2 "$SC/faser2_senscalc_fpf.dat" ) > "$W/bc10_published_table.log" 2>&1 \
    || { echo "FAIL Fig. 7 (published table, rebuilt curves)"; rc=1; }
fi

# --- results/: the rebuilt tables, the rebuilt SensCalc curves and the validation reports
put() {  # put <src> <dst>: copy if the source exists
  if [ -s "$1" ]; then mkdir -p "$(dirname "$2")"; cp -p "$1" "$2"; else echo "MISSING $1" >&2; rc=1; fi
}
for b in bc4:work/scalar/results/bc4 bc5:work/scalar/results/bc5 bc10:work/alp/tables_main/bc10 hnl:work/hnl/table/hnl; do
  for f in sensitivity.csv sensitivity_nsig10.csv; do put "${b#*:}/$f" "results/tables/${b%%:*}/$f"; done
done
for f in ship_senscalc.dat faser2_senscalc_fpf.dat; do put "$SC/$f" "results/senscalc/$f"; done
SANITIZE=(-e "s#$REPRO/#./#g" -e "s#$REPRO#.#g")
case "${GRENDEL_REF:-none}" in none) ;; *) SANITIZE+=(-e "s#$GRENDEL_REF#\$GRENDEL_REF#g") ;; esac
SANITIZE+=(-e "s#$HOME#~#g")
for p in fonll scalar hnl alp senscalc; do
  [ -d "work/$p/validation" ] || continue
  rm -rf "results/validation/$p"; mkdir -p "results/validation/$p"
  (cd "work/$p/validation" && find . -type f) | while read -r f; do
    mkdir -p "results/validation/$p/$(dirname "$f")"
    if grep -qI . "work/$p/validation/$f"; then   # text: paths of this machine -> ".", $GRENDEL_REF, ~
      sed "${SANITIZE[@]}" "work/$p/validation/$f" > "results/validation/$p/$f"
    else
      cp -p "work/$p/validation/$f" "results/validation/$p/$f"
    fi
  done
done

# --- comparison with the paper's PNGs
{
  echo "Rebuilt Sec. IV figures vs the paper's PNGs ($PAPER), $(date '+%F %T')"
  echo "A pixel differs when any channel differs; '>32/255' counts the visible differences."
  cmp_one() {  # cmp_one <label> <ours.png> <paper name>
    echo; echo "== $1"
    if [ -s "$2" ]; then python "$CMP" "$2" "$PAPER/$3" "work/figures/diff/$(basename "$(dirname "$2")")_$3"
    else echo "not built: $2"; fi
  }
  cmp_one "Fig. 4 (BC4)"  figures/scalar/bc4_grendel_paper.png      bc4_grendel_paper.png
  cmp_one "Fig. 5 (BC5)"  figures/scalar/bc5_grendel_paper.png      bc5_grendel_paper.png
  cmp_one "Fig. 6 (HNL)"  figures/hnl/hnlimits_grendel_paper.png    hnlimits_grendel_paper.png
  cmp_one "Fig. 7 (BC10), rebuilt table and rebuilt SensCalc curves" figures/alp/bc10_grendel_paper.png bc10_grendel_paper.png
  cmp_one "Fig. 7 (BC10), rebuilt table, tracked SensCalc curves"    work/figures/bc10_tracked_curves/bc10_grendel_paper.png bc10_grendel_paper.png
  cmp_one "Fig. 7 (BC10), published table, rebuilt SensCalc curves"  work/figures/bc10_published_table/bc10_grendel_paper.png bc10_grendel_paper.png
  echo
  echo "== Input digests (sha256, first 16 hex digits)"
  for f in figures/scalar/*.png figures/hnl/*.png figures/alp/*.png results/senscalc/ship_senscalc.dat \
           results/senscalc/faser2_senscalc_fpf.dat \
           "$SRC/curves/data/bc10/projections/ship_senscalc.dat" "$SRC/curves/data/bc10/projections/faser2_senscalc_fpf.dat"; do
    [ -s "$f" ] && printf '%s  %s\n' "$(sha256sum "$f" | cut -c1-16)" "${f#"$REPRO/"}"
  done
} > "$REPORT"
cat "$REPORT"

need figures/scalar/bc4_grendel_paper.png "pipelines/scalar step 10"
need figures/scalar/bc5_grendel_paper.png "pipelines/scalar step 10"
need figures/hnl/hnlimits_grendel_paper.png "pipelines/hnl step 11"
need figures/alp/bc10_grendel_paper.png "figures/collect.sh"
exit $rc
