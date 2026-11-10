#!/bin/bash
# HNL step 11: Fig. 6 (hnlimits_grendel_paper.png) from a results tree, with the paper's renderer:
#     python -m curves.plot paper --grendel-dir <dir> --benchmark hnl      (needs HNLimits 1.2.0)
# Renders both the rebuilt table (work/hnl/table) and the published one (curves/data/grendel_v2_3) into
# work/hnl/figures/{rebuilt,published}/, copies the rebuilt PNG and PDF to figures/hnl/ and compares both PNGs, pixel
# by pixel, with the paper's figure ($GRENDEL_PAPER_FIGURES = reference/paper_figures, in this repository).
#     11_plot.sh [rebuilt|published|both]
set -euo pipefail
source "$(dirname "$0")/../../env/activate.sh"
F=$REPRO/work/hnl/figures; mkdir -p $F $REPRO/figures/hnl
export MPLCONFIGDIR=$F/.mpl PYTHONDONTWRITEBYTECODE=1
cd $REPRO/src/grendel-sensitivity
plot() {  # plot <name> <results tree with hnl/sensitivity.csv>
  python $REPRO/pipelines/common/slot_run.py -n 1 --label hnl:11_plot_$1 -- \
    python -m curves.plot paper --grendel-dir $2 --benchmark hnl --out $F/$1 > $REPRO/logs/hnl/11_plot_$1.log 2>&1
  ls $F/$1
}
what=${1:-both}
if [ $what != rebuilt ]; then
  P=$F/published_tree; mkdir -p $P/hnl
  cp $GRENDEL_PUBLISHED/grendel_v2_3/hnl/sensitivity*.csv $P/hnl/
  plot published $P
fi
if [ $what != published ]; then
  [ -s $REPRO/work/hnl/table/hnl/sensitivity.csv ] || { echo "no rebuilt table yet (step 09)"; exit 1; }
  plot rebuilt $REPRO/work/hnl/table
  cp $F/rebuilt/hnlimits_grendel_paper.png $REPRO/figures/hnl/hnlimits_grendel_paper.png
  cp $F/rebuilt/hnlimits_grendel_paper.pdf $REPRO/figures/hnl/hnlimits_grendel_paper.pdf
fi
# pixel comparison with the paper's figure (the Overleaf PNG)
O=$GRENDEL_PAPER_FIGURES/hnlimits_grendel_paper.png
V=$REPRO/work/hnl/validation; mkdir -p $V
if [ -f $O ]; then
  for n in published rebuilt; do
    [ -f $F/$n/hnlimits_grendel_paper.png ] || continue
    python $REPRO/pipelines/hnl/lib/compare_png.py $F/$n/hnlimits_grendel_paper.png $O --diff $V/figure_${n}_diff.png \
      > $V/figure_$n.txt
    cat $V/figure_$n.txt
  done
else
  echo "missing $O (reference/paper_figures is part of the repository): no pixel comparison"
fi
