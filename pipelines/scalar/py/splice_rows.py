"""Our tables with the rows of some masses replaced by the published lines (a figure-attribution check of step 10).

    python splice_rows.py OURS_DIR PUBLISHED_DIR OUT_DIR MASS [MASS ...]

For bc4 and bc5, sensitivity.csv and sensitivity_nsig10.csv: every text line is ours, except the lines of the given
masses, which are copied from the published table. Rendered, the result shows whether the remaining rows (ours)
produce the paper's pixels: 10_figures.sh passes the masses whose decay templates differ from the published ones.
"""
import sys
from pathlib import Path

ours, pub, out = (Path(p) for p in sys.argv[1:4])
split = {float(m) for m in sys.argv[4:]}
for bench in ("bc4", "bc5"):
    (out / bench).mkdir(parents=True, exist_ok=True)
    for table in ("sensitivity.csv", "sensitivity_nsig10.csv"):
        lo = (ours / bench / table).read_text().splitlines()
        lp = (pub / bench / table).read_text().splitlines()
        assert lo[0] == lp[0], f"{bench}/{table}: headers differ"
        i = lo[0].split(",").index("mass_GeV")
        mass = lambda line: float(line.split(",")[i])
        published = {mass(line): line for line in lp[1:]}
        lines = [published[mass(line)] if mass(line) in split else line for line in lo[1:]]
        (out / bench / table).write_text("\n".join([lo[0]] + lines) + "\n")
        print(f"{bench}/{table}: {sum(mass(line) in split for line in lo[1:])} rows from the published table")
