#!/usr/bin/env bash
# Step 09 -- the tables of curves/data/grendel_v2_5/{bc4,bc5} from the one-mass scans, the Sec. IV numbers, and
# the row-by-row comparison with the published tables.
#
# Published recipe: scalar_eventcalc_20261019/collect.py (copied to py/ref/): concatenate the per-job rows sorted by
# mass, `python -m grendel.io.thresholds sensitivity.csv --threshold 10`, one BC5 run_metadata.json; BC5
# paper_numbers.txt = bc5_v2_1_20261003/results/bc5_hybrid/paper_numbers.py on the table (py/ref/paper_numbers.py,
# unchanged apart from its default CURVES_ROOT); BC4 Sec. IV numbers = v2_3_20261016/edge_ratios.py (py/ref/edge_ratios.py, unchanged).
# Output: work/scalar/results/{bc4,bc5}/ (the layout `curves.plot paper --grendel-dir` reads) and
# work/scalar/validation/{final_bc4,final_bc5,sec4_numbers}.txt. Seconds. Runs on partial results too (then it
# says which masses are missing and exits 1 after writing what it has).
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
step_lock
R="$W/results"
mkdir -p "$R" "$W/validation"
rc=0
# masses whose decay templates are not event by event the published ones (06's report; their rows are compared apart).
# 06 writes it only with the reference data ($GRENDEL_REF); without it every mass is compared in one group.
SPLIT=""
if [ -s "$W/validation/templates.tsv" ]; then
  SPLIT=$(awk -F'\t' 'NR > 1 && $5 != 1 {print $1}' "$W/validation/templates.tsv" | sort -g | paste -sd ' ')
else
  echo "NOTE no template report from 06 (reference data not available): the rows are not split by template identity" >&2
fi
python "$S/py/collect.py" "$W" main "$R" "$BC4_MASSES" "$BC5_MASSES" "$SRC" || rc=1
for b in bc4 bc5; do
  [ -s "$R/$b/sensitivity.csv" ] || continue
  # per row (ratios of the edges, flags) and per column (bit-identical columns, largest relative difference)
  python "$S/py/compare_rows.py" "$b" "$PUB/$b/sensitivity.csv" "$R/$b/sensitivity.csv" \
    --tsv "$W/validation/final_${b}.tsv" > "$W/validation/final_${b}.txt"
  tail -n 4 "$W/validation/final_${b}.txt"
  # shellcheck disable=SC2086
  python "$S/py/compare_tables.py" "$PUB/$b/sensitivity.csv" "$R/$b/sensitivity.csv" \
    --tsv "$W/validation/final_${b}_columns.tsv" --split $SPLIT > "$W/validation/final_${b}_columns.txt"
  grep -E "SHA-256|bit-identical|largest|verdict|has_sensitivity:|has_sensitivity_N10:" \
    "$W/validation/final_${b}_columns.txt" | cut -c1-200
  # sensitivity_nsig10.csv is grendel.io.thresholds' projection of sensitivity.csv; it re-parses the table with
  # pandas' default float parser, which differs between the published machine (pandas 3.0.6, macOS arm64) and the
  # venv here (3.1.1: correctly rounded) in the last 1-3 of 17 digits. So ours is compared with the published file
  # and with the published sensitivity.csv projected here (the latter isolates the data from the parser).
  python -m grendel.io.thresholds "$PUB/$b/sensitivity.csv" --threshold 10 \
    --out "$W/validation/published_${b}_nsig10_projected_here.csv" > /dev/null
  # shellcheck disable=SC2086
  { python "$S/py/compare_tables.py" "$PUB/$b/sensitivity_nsig10.csv" "$R/$b/sensitivity_nsig10.csv" --split $SPLIT
    python "$S/py/compare_tables.py" "$W/validation/published_${b}_nsig10_projected_here.csv" \
      "$R/$b/sensitivity_nsig10.csv" --split $SPLIT; } > "$W/validation/final_${b}_nsig10.txt"
done
# BC5 run_metadata.json: every field against the published one (the input digests differ by construction, README)
[ -s "$R/bc5/run_metadata.json" ] && python - "$PUB/bc5/run_metadata.json" "$R/bc5/run_metadata.json" <<'EOF' \
  > "$W/validation/final_bc5_run_metadata.txt"
import json, sys
p, o = (json.load(open(f)) for f in sys.argv[1:3])
def flat(d, pre=""):
    for k, v in d.items():
        yield from (flat(v, f"{pre}{k}.") if isinstance(v, dict) and k != "inputs" else [(pre + k, v)])
P, O = dict(flat(p)), dict(flat(o))
same = [k for k in P if k in O and P[k] == O[k]]
print(f"run_metadata.json: {len(same)}/{len(P)} fields equal: {' '.join(same)}")
for k in sorted(set(P) | set(O)):
    if k == "inputs":
        a, b = P.get(k, {}), O.get(k, {})
        print(f"  inputs: {len(b)} per-mass digests here, {len(a)} published, {sum(a.get(m) == b[m] for m in b)} equal")
    elif P.get(k, "<absent>") != O.get(k, "<absent>"):
        s = lambda x: (json.dumps(x)[:110] + "...") if len(json.dumps(x)) > 110 else json.dumps(x)
        print(f"  {k}: published {s(P.get(k, '<absent>'))}, here {s(O.get(k, '<absent>'))}")
EOF
{
  if [ -s "$R/bc5/sensitivity.csv" ]; then
    if (cd "$R/.." && python "$S/py/ref/paper_numbers.py" results/bc5/sensitivity.csv "$SRC") > "$R/bc5/paper_numbers.txt" 2>/dev/null; then
      echo "== BC5 paper_numbers.txt: ours vs published (first line is the table path)"
      diff "$R/bc5/paper_numbers.txt" "$PUB/bc5/paper_numbers.txt" && echo "identical" || true
    else
      rm -f "$R/bc5/paper_numbers.txt"; echo "== BC5 paper_numbers.py needs the complete table"
    fi
  fi
  if [ -s "$R/bc4/sensitivity.csv" ]; then
    # (full outputs to files: a `| head` here would SIGPIPE the script under pipefail)
    PYTHONPATH="$SRC" python "$S/py/ref/edge_ratios.py" bc4 "$R/bc4/sensitivity.csv" \
      > "$W/validation/bc4_edge_ratios_ours.txt" 2>/dev/null || true
    PYTHONPATH="$SRC" python "$S/py/ref/edge_ratios.py" bc4 "$PUB/bc4/sensitivity.csv" \
      > "$W/validation/bc4_edge_ratios_published.txt" 2>/dev/null || true
    echo "== BC4 Sec. IV numbers (edge_ratios.py): ours vs published (validation/bc4_edge_ratios_*.txt)"
    diff "$W/validation/bc4_edge_ratios_ours.txt" "$W/validation/bc4_edge_ratios_published.txt" && echo "identical" || true
  fi
  echo "== paper Sec. IV: deepest BC4 7.2e-12 at 1.0 GeV, BC4 reach 3.6 GeV, BC5 2.0e-17 near 50 GeV"
  python - "$R" <<'EOF'
import sys
from pathlib import Path
import pandas as pd
R = Path(sys.argv[1])
for b in ("bc4", "bc5"):
    p = R / b / "sensitivity.csv"
    if not p.exists():
        continue
    t = pd.read_csv(p)
    s = t[t.has_sensitivity.astype(bool)]
    print(f"{b}: deepest {s.u2_min.min():.2g} at {s.loc[s.u2_min.idxmin(), 'mass_GeV']:g} GeV; reach (N>=3) "
          f"{s.mass_GeV.max():g} GeV; {len(s)} sensitive masses of {len(t)}; N>=10 {int(t.has_sensitivity_N10.sum())}")
EOF
} 2>&1 | tee "$W/validation/sec4_numbers.txt"
[ $rc -eq 0 ] && echo "09_collect OK" || echo "09_collect: tables incomplete"
exit $rc
