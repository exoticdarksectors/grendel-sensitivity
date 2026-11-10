"""Smoke run of the whole central-grid chain on a tiny grid (minutes instead of hours).

    source env/activate.sh && source pipelines/fonll/env.sh
    python pipelines/common/slot_run.py -n 1 --label fonll:smoke -- python pipelines/fonll/smoke.py --quark charm

FONLL has no statistics knob (its integration accuracy is internal), so "tiny statistics" means a tiny grid: the
module-level grid of grendel.production.fonll_grids.generate is shrunk in place before the unchanged generate()
runs, to 4 rapidities x 5 meson pT nodes on a 20-node raw quark pT grid (instead of 100 x 100 on 80). Everything
else is the production path: fonllgridlha (rapidity chunks), the charm CTEQ6.6 feeddown calibration and its fit,
both fragmfonll passes, the D* -> D0 combination and the grid writer. Outputs go to work/fonll/smoke/ and are
not physics results.
"""
import argparse
import sys
from pathlib import Path

from grendel.production.fonll_grids import generate as gen

ap = argparse.ArgumentParser()
ap.add_argument("--quark", action="append", choices=sorted(gen.QUARKS))
ap.add_argument("--grid-workers", type=int, default=1)
ap.add_argument("--root", type=Path, required=True)
args = ap.parse_args()

# shrink the grid in place: Y_VALUES/PT_VALUES are the list objects every function (and grid_input's default)
# refers to; the raw-grid pT node count is an int default, so grid_input's defaults are rebound as well.
gen.Y_VALUES[:] = [-3.0, -1.0, 1.0, 3.0]
gen.PT_VALUES[:] = [0.0, 2.0, 5.0, 20.0, 50.0]
gen.FONLL_GRID_NPT = 20
gen.grid_input.__defaults__ = (1.0, 1.0, gen.Y_VALUES, gen.FONLL_GRID_PTMAX, gen.FONLL_GRID_NPT)

ws = gen.Workspace.resolve(None, args.root / "output", args.root / "run", args.root / "logs")
ws.out.mkdir(parents=True, exist_ok=True)
ok = True
for quark in args.quark or ["bottom", "charm"]:
    calib = gen.load_or_fit_charm_feeddown_weight(ws) if quark == "charm" else None
    if calib:
        print(f"[smoke charm] feeddown weight {calib['vector_to_direct_weight']:.10g} "
              f"(max |rel| vs public {calib['max_abs_relative_difference']:.2e}; production fit: 1.21525336412)")
    s = gen.generate(ws, "nlo", quark, gen.central_variation(quark), reuse_existing_grids=False,
                     grid_workers=args.grid_workers, feeddown_calibration=calib)
    rows_ok = s["rows"] == len(gen.Y_VALUES) * len(gen.PT_VALUES)
    print(f"[smoke {quark}] {s['path']}: rows={s['rows']} min={s['min_dsigma_dpT_dy_pb_per_GeV']:.4e} "
          f"max={s['max_dsigma_dpT_dy_pb_per_GeV']:.4e} trapz={s['trapezoid_integral_pb']:.4e} pb "
          f"{'OK' if rows_ok and s['max_dsigma_dpT_dy_pb_per_GeV'] > 0 else 'FAIL'}")
    ok = ok and rows_ok
sys.exit(0 if ok else 1)
