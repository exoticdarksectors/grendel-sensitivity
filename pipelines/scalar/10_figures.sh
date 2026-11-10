#!/usr/bin/env bash
# Step 10 -- Figs. 4 and 5: render our tables and the published ones with the same renderer and compare both
# pixel by pixel with the paper's PNGs (reference/paper_figures/bc{4,5}_grendel_paper.png, from the Overleaf project).
#
# Published recipe: python -m curves.plot paper --grendel-dir curves/data/grendel_v2_5 --benchmark bc4 bc5
# (curves/data/grendel_v2_5/README.md). The renderer is that of 8802f50 (the paper figures are its output).
# Attribution check: our tables with the rows of the masses whose decay templates differ from the published ones
# (06's report) replaced by the published rows (py/splice_rows.py), rendered the same way (figures/spliced): if that
# is pixel-identical to the paper, every pixel difference of ours comes from those few template rows. It needs 06's
# template report, i.e. the reference data ($GRENDEL_REF), and is skipped without it. The PDFs are also compared with
# the paper's PDFs (reference/paper_figures).
# Output: work/scalar/figures/{ours,published_tables,spliced}/bc{4,5}_grendel_paper.{png,pdf}, diff images, and
# work/scalar/validation/figures.txt; our figures are copied to figures/scalar/. Seconds.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
step_lock
F="$W/figures"
PAPER="$GRENDEL_PAPER_FIGURES"                       # the paper's PNGs (in this repository)
PAPER_PDF="$GRENDEL_PAPER_FIGURES"                   # the paper's PDFs (in this repository)
mkdir -p "$F/ours" "$F/published_tables" "$W/validation"
cd "$SRC"
python -m curves.plot paper --grendel-dir "$PUB" --benchmark bc4 bc5 --out "$F/published_tables" >> "$L/10_figures.log" 2>&1
python -m curves.plot paper --grendel-dir "$W/results" --benchmark bc4 bc5 --out "$F/ours" >> "$L/10_figures.log" 2>&1
SPLICE=0
if ref_ok "spliced figures (needs 06's template report)" "$W/validation/templates.tsv"; then
  SPLICE=1
  SPLIT=$(awk -F'\t' 'NR > 1 && $5 != 1 {print $1}' "$W/validation/templates.tsv" | sort -g | paste -sd ' ')
  # shellcheck disable=SC2086
  python "$S/py/splice_rows.py" "$W/results" "$PUB" "$F/spliced_tables" $SPLIT >> "$L/10_figures.log"
  python -m curves.plot paper --grendel-dir "$F/spliced_tables" --benchmark bc4 bc5 --out "$F/spliced" \
    >> "$L/10_figures.log" 2>&1
fi
PDFS=0
ref_ok "PDFs vs the Overleaf PDFs" "$PAPER_PDF/bc4_grendel_paper.pdf" "$PAPER_PDF/bc5_grendel_paper.pdf" && PDFS=1
{
  for b in bc4 bc5; do
    echo "== $b: published tables rendered here vs the paper PNG (renderer check)"
    python "$S/py/compare_png.py" "$F/published_tables/${b}_grendel_paper.png" "$PAPER/${b}_grendel_paper.png" \
      "$F/published_tables/${b}_diff_vs_paper.png"
    echo "== $b: our tables vs the paper PNG (data check)"
    python "$S/py/compare_png.py" "$F/ours/${b}_grendel_paper.png" "$PAPER/${b}_grendel_paper.png" \
      "$F/ours/${b}_diff_vs_paper.png"
    if [ "$SPLICE" -eq 1 ]; then
      echo "== $b: our tables with the published rows at the template-differing masses ($SPLIT GeV) vs the paper PNG"
      python "$S/py/compare_png.py" "$F/spliced/${b}_grendel_paper.png" "$PAPER/${b}_grendel_paper.png"
    fi
    [ "$PDFS" -eq 1 ] || continue
    echo "== $b: PDFs vs the paper PDF (bytes, apart from the /CreationDate stamp)"
    for d in published_tables ours; do
      python - "$F/$d/${b}_grendel_paper.pdf" "$PAPER_PDF/${b}_grendel_paper.pdf" "$d" <<'EOF'
import re, sys
a, b = (re.sub(rb"/CreationDate \(D:[^)]*\)", b"", open(p, "rb").read()) for p in sys.argv[1:3])
print(f"{sys.argv[3]}: {'identical' if a == b else 'DIFFERENT'} ({len(a)} / {len(b)} bytes without the date stamp)")
EOF
    done
  done
} | tee "$W/validation/figures.txt"
mkdir -p "$REPRO/figures/scalar"
cp "$F"/ours/bc[45]_grendel_paper.png "$F"/ours/bc[45]_grendel_paper.pdf "$REPRO/figures/scalar/"
step_ok 10_figures
