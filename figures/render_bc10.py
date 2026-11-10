"""Render Fig. 7 (BC10) from the rebuilt table with the rebuilt SensCalc curves.

    python figures/render_bc10.py TABLES_DIR OUT_DIR [--ship FILE] [--faser2 FILE]

``curves.plot paper --benchmark bc10`` draws the SHiP and FASER2 projections from the tracked files
``src/grendel-sensitivity/curves/data/bc10/projections/{ship_senscalc,faser2_senscalc_fpf}.dat``.
To draw the figure entirely from rebuilt inputs, this script points those two entries of
``curves.plot.bc_panel.BENCHMARKS["bc10"]`` at the files written by ``pipelines/senscalc``
(``work/senscalc/out/projections/``) and then runs the unchanged renderer. Nothing in
``src/`` is modified. Every other curve (existing bounds, DarkQuest, CODEX-b) is a digitised
literature curve and is drawn from the tracked file, as in the paper.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tables", type=Path, help="results tree holding bc10/sensitivity.csv (and _nsig10)")
    ap.add_argument("out", type=Path)
    ap.add_argument("--ship", type=Path, help="replacement for projections/ship_senscalc.dat")
    ap.add_argument("--faser2", type=Path, help="replacement for projections/faser2_senscalc_fpf.dat")
    a = ap.parse_args()

    from curves.plot import bc_panel
    from curves.plot.__main__ import main as plot_main

    swap = {"ship_senscalc.dat": a.ship, "faser2_senscalc_fpf.dat": a.faser2}
    done = set()
    for entry in bc_panel.BENCHMARKS["bc10"]["projections"]:
        new = swap.get(Path(entry["path"]).name)
        if new is not None:
            if not new.is_file():
                print(f"missing {new}", file=sys.stderr)
                return 2
            print(f"projection {Path(entry['path']).name} <- {new}")
            entry["path"] = new.resolve()
            done.add(Path(new).name)
    wanted = {p.name for p in swap.values() if p is not None}
    if wanted - done:
        print(f"no projection entry matched {sorted(wanted - done)}", file=sys.stderr)
        return 2
    return plot_main(["paper", "--grendel-dir", str(a.tables), "--benchmark", "bc10", "--out", str(a.out)])


if __name__ == "__main__":
    sys.exit(main())
