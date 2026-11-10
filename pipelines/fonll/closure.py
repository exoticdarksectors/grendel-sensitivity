"""Bottom public-point closure: CTEQ6.6 B hadron at pT = 5 GeV, y = 0 (step 02b).

Same code path as the NNPDF4.0 central bottom grid, only the LHAPDF set differs (10550 instead of 331700):
fonllgridlha at y = -0.5, 0, 0.5 (80 raw pT nodes up to 200 GeV, the production card), then fragmfonll with
Kartvelishvili alpha = 24.2 at (5, 0). Compared with the public FONLL v1.3.2 web form (7.5329e6 pb/GeV,
grendel.production.fonll_grids.validate.PUBLIC_BOTTOM_PT5_Y0) and with the original macOS run (7.535277e6,
data/fonll/PROVENANCE.md).
"""
import argparse
import json
from pathlib import Path

from grendel.production.fonll_grids import generate as gen
from grendel.production.fonll_grids.validate import PUBLIC_BOTTOM_PT5_Y0

ORIGINAL_LOCAL = 7.535277e6

ap = argparse.ArgumentParser()
ap.add_argument("--out", type=Path, required=True)
args = ap.parse_args()

ws = gen.Workspace.resolve()          # GRENDEL_FONLL_DIR / GRENDEL_FONLL_GRID_DIR from pipelines/fonll/env.sh
run_dir = ws.run / "cteq66_bottom_y0"
run_dir.mkdir(parents=True, exist_ok=True)
grid = run_dir / "ct_b_.out"
ys = [-0.5, 0.0, 0.5]
if not (grid.exists() and len(grid.read_text().splitlines()) == 1 + len(ys) * gen.FONLL_GRID_NPT):
    for f in run_dir.glob("ct_b_*"):
        f.unlink()
    gen.run_command([str(ws.executable("fonllgridlha"))], run_dir, gen.grid_input("ct_b_", 10550, 4.75, y_values=ys),
                    ws.logs / "closure_cteq66_bottom_y0_fonllgrid.log", ws.env())
b = gen.QUARKS["bottom"]
values, _ = gen.run_fragmentation(ws, "cteq66", "bottom", run_dir, grid, int(b["frag_mode"]), float(b["frag_param"]),
                                  [(5.0, 0.0)], "public_point_pt5_y0")
local = values[gen.point_key(5.0, 0.0)]
res = {"point": "B hadron, pT = 5 GeV, y = 0, CTEQ6.6, Kartvelishvili alpha = 24.2",
       "local_linux": local, "public_fonll_web": PUBLIC_BOTTOM_PT5_Y0, "original_local_macos": ORIGINAL_LOCAL,
       "rel_vs_public": (local - PUBLIC_BOTTOM_PT5_Y0) / PUBLIC_BOTTOM_PT5_Y0,
       "rel_vs_original_local": (local - ORIGINAL_LOCAL) / ORIGINAL_LOCAL, "raw_grid": str(grid)}
args.out.parent.mkdir(parents=True, exist_ok=True)
args.out.write_text(json.dumps(res, indent=2) + "\n")
print(json.dumps(res, indent=2))
