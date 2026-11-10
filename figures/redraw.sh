#!/usr/bin/env bash
# Redraw the four Sec. IV figures from stored tables (seconds) and compare them pixel by pixel with the reference
# PNGs. Needs only the analysis venv: bash env/bootstrap.sh --minimal && source env/activate.sh.
#
#   bash figures/redraw.sh published   the published tables (src/grendel-sensitivity/curves/data)
#                                      -> the paper's figures (reference/paper_figures)
#   bash figures/redraw.sh rebuilt     the tables rebuilt by the pipelines (results/tables, results/senscalc)
#                                      -> the rebuilt figures (figures/)
#
# Output in work/redraw/<published|rebuilt>/; exit status 0 if every figure is pixel-identical to its reference.
set -uo pipefail
REPRO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [ -z "${GRENDEL_PUBLISHED:-}" ]; then
  # shellcheck disable=SC1091
  source "$REPRO/env/activate.sh"
fi
WHICH="${1:-published}"
SRC="$REPRO/src/grendel-sensitivity"
CMP="$REPRO/pipelines/scalar/py/compare_png.py"
OUT="$REPRO/work/redraw/$WHICH"
export MPLCONFIGDIR="$REPRO/work/redraw/.mpl" PYTHONDONTWRITEBYTECODE=1
mkdir -p "$OUT" "$MPLCONFIGDIR"

render() {  # render <benchmark> <tables dir> [--ship F --faser2 F]
  local b=$1 t=$2
  shift 2
  if [ $# -gt 0 ]; then
    (cd "$SRC" && python "$REPRO/figures/render_bc10.py" "$t" "$OUT/$b" "$@")
  else
    (cd "$SRC" && python -m curves.plot paper --grendel-dir "$t" --benchmark "$b" --out "$OUT/$b")
  fi > "$OUT/$b.log" 2>&1 || { echo "FAIL rendering $b (log: ${OUT#"$REPRO"/}/$b.log)"; return 1; }
}

case "$WHICH" in
  published)
    P="$GRENDEL_PUBLISHED"
    render bc4 "$P/grendel_v2_5" && render bc5 "$P/grendel_v2_5" && render hnl "$P/grendel_v2_3" \
      && render bc10 "$P/grendel_v2_4" || exit 1
    ref() { echo "$GRENDEL_PAPER_FIGURES/$1"; } ;;
  rebuilt)
    T="$REPRO/results/tables"
    render bc4 "$T" && render bc5 "$T" && render hnl "$T" \
      && render bc10 "$T" --ship "$REPRO/results/senscalc/ship_senscalc.dat" \
           --faser2 "$REPRO/results/senscalc/faser2_senscalc_fpf.dat" || exit 1
    ref() { case "$1" in bc10*) echo "$REPRO/figures/alp/$1" ;; hnl*) echo "$REPRO/figures/hnl/$1" ;;
                        *) echo "$REPRO/figures/scalar/$1" ;; esac; } ;;
  *) echo "usage: redraw.sh published|rebuilt" >&2; exit 2 ;;
esac

rc=0
cd "$REPRO"
for f in bc4/bc4_grendel_paper.png bc5/bc5_grendel_paper.png hnl/hnlimits_grendel_paper.png bc10/bc10_grendel_paper.png; do
  r="$(ref "$(basename "$f")")"
  if cmp -s "$OUT/$f" "$r"; then
    echo "identical (bytes)   ${OUT#"$REPRO"/}/$f = ${r#"$REPRO"/}"
  else
    line=$(python "$CMP" "$OUT/$f" "$r" | tail -1)
    case "$line" in
      "pixels differing: 0.00000%"*) echo "identical (pixels)  ${OUT#"$REPRO"/}/$f = ${r#"$REPRO"/}" ;;
      *) echo "DIFFERENT           ${OUT#"$REPRO"/}/$f vs ${r#"$REPRO"/}: $line"; rc=1 ;;
    esac
  fi
done
exit $rc
