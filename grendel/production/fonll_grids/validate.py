"""Closure checks for the grid generator.

    python -m grendel.production.fonll_grids.validate gate [--quark bottom] [--rtol 1e-3] [--grid-workers 6]
    python -m grendel.production.fonll_grids.validate public-points [--tolerance 0.05]
    python -m grendel.production.fonll_grids.validate charm-point

``gate`` regenerates the central grid (member 0, scale (1,1), central
mass) and asserts its dsigma column reproduces the tracked reference in
``grendel/production/data/fonll/central`` within ``--rtol``; run it before
a long campaign so days of compute are not spent on a moved central. It
also reports the fraction of the cross section in the outermost pT and
rapidity bins as a grid-coverage check.

``public-points`` compares a CTEQ6.6 bottom grid produced by this code
(``generate --pdf cteq66 --quark bottom``) at pT = 5 GeV, y = 0 with the
public FONLL v1.3.2 web form, and refits the charm D* feeddown weight
against the nine public D0 points. ``charm-point`` is the standalone
direct-BCFY pseudoscalar D0 point at pT = 5 GeV, y = 0 (no feeddown), the
+4.5 % diagnostic the feeddown calibration absorbs.

Exit status is nonzero when a check fails.
"""
from __future__ import annotations

import argparse
import math
import sys
import tempfile
from pathlib import Path

import numpy as np

from . import generate as gen

REFERENCE_DIR = Path(__file__).resolve().parents[1] / "data" / "fonll" / "central"

# Public FONLL web v1.3.2, CTEQ6.6, pp 14 TeV, central, dsigma/dpT/dy at
# pT = 5 GeV, y = 0, fragmentation fraction 1; queried 2026-05-29.
PUBLIC_BOTTOM_PT5_Y0 = 7.5329e6      # B hadron default, N=5 Kartvelishvili alpha=24.2
PUBLIC_CHARM_D0_PT5_Y0 = 3.5221e7    # D0 default, BCFY r=0.1 including the D* feeddown convention


def data_columns(path: Path):
    data = np.loadtxt(path, comments="#")
    return data[:, 0], data[:, 1], data[:, 2]


def coverage_tail_bound(ds: np.ndarray) -> dict:
    n_pt, n_y = len(gen.PT_VALUES), len(gen.Y_VALUES)
    grid = ds.reshape(n_pt, n_y)
    total = float(np.trapezoid(np.trapezoid(grid, gen.Y_VALUES, axis=1), gen.PT_VALUES))
    pt_tail = float(np.trapezoid(grid[-1, :], gen.Y_VALUES)) * (gen.PT_VALUES[-1] - gen.PT_VALUES[-2])
    y_tail = (float(np.trapezoid(grid[:, 0], gen.PT_VALUES)) + float(np.trapezoid(grid[:, -1], gen.PT_VALUES))) \
        * (gen.Y_VALUES[1] - gen.Y_VALUES[0])
    return {"total_pb": total, "pt_edge_fraction": pt_tail / total if total else float("nan"),
            "y_edge_fraction": y_tail / total if total else float("nan")}


def gate_quark(ws: gen.Workspace, quark: str, reference_dir: Path, rtol: float, grid_workers: int) -> bool:
    feeddown = gen.load_or_fit_charm_feeddown_weight(ws) if quark == "charm" else None
    summary = gen.generate(ws, "nlo", quark, gen.central_variation(quark), reuse_existing_grids=False,
                           grid_workers=grid_workers, feeddown_calibration=feeddown)
    _, _, fresh = data_columns(Path(summary["path"]))
    _, _, ref = data_columns(reference_dir / gen.grid_filename("nlo", "central", quark))
    if fresh.shape != ref.shape:
        print(f"[{quark}] FAIL: shape {fresh.shape} vs reference {ref.shape}")
        return False
    nz = ref != 0
    rel = np.zeros_like(ref)
    rel[nz] = np.abs(fresh[nz] - ref[nz]) / np.abs(ref[nz])
    max_rel = float(rel.max())
    cov = coverage_tail_bound(fresh)
    print(f"[{quark}] {'PASS' if max_rel <= rtol else 'FAIL'}: max relative diff vs tracked central = "
          f"{max_rel:.3e} (rtol={rtol:.1e})")
    print(f"[{quark}] coverage: total={cov['total_pb']:.6e} pb; pT-edge fraction={cov['pt_edge_fraction']:.3e}; "
          f"y-edge fraction={cov['y_edge_fraction']:.3e}")
    return max_rel <= rtol


def find_point(path: Path, pt: float, y: float) -> float:
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.split()
        if len(cols) < 3:
            continue
        if (math.isclose(float(cols[0].replace("D", "E")), pt, abs_tol=1e-12)
                and math.isclose(float(cols[1].replace("D", "E")), y, abs_tol=1e-12)):
            return float(cols[2].replace("D", "E"))
    raise ValueError(f"point pT={pt}, y={y} not found in {path}")


def public_points(ws: gen.Workspace, tolerance: float) -> bool:
    ok = True
    path = ws.out / gen.grid_filename("cteq66", "central", "bottom")
    if not path.exists():
        raise FileNotFoundError(f"{path} not found; run generate --pdf cteq66 --quark bottom first")
    local = find_point(path, 5.0, 0.0)
    rel = (local - PUBLIC_BOTTOM_PT5_Y0) / PUBLIC_BOTTOM_PT5_Y0
    status = "OK" if abs(rel) <= tolerance else "FAIL"
    print(f"{status} bottom: local={local:.6e} public={PUBLIC_BOTTOM_PT5_Y0:.6e} rel_diff={rel:+.3%}")
    ok = ok and status == "OK"

    charm = gen.fit_charm_d0_feeddown_weight(ws)
    max_abs_rel = float(charm["max_abs_relative_difference"])
    status = "OK" if max_abs_rel <= tolerance else "FAIL"
    print(f"{status} charm: calibrated public-D0 feeddown max_abs_rel_diff={max_abs_rel:.3%} over "
          f"{len(charm['points'])} CTEQ6.6 y=0 public points (weight {charm['vector_to_direct_weight']:.10g})")
    for point in charm["points"]:
        if math.isclose(point["pt"], 5.0, abs_tol=1e-12):
            print(f"  pT=5,y=0 local={point['local']:.6e} public={point['public']:.6e} "
                  f"rel_diff={point['relative_difference']:+.3%}")
    return ok and status == "OK"


def charm_point(ws: gen.Workspace) -> None:
    grid_file = gen.ensure_charm_calibration_grid(ws)
    values, _ = gen.run_fragmentation(ws, "cteq66", "charm", grid_file.parent, grid_file, 5, 0.1, [(5.0, 0.0)],
                                      "direct_pseudoscalar_pt5_y0")
    local = values[gen.point_key(5.0, 0.0)]
    rel = (local - PUBLIC_CHARM_D0_PT5_Y0) / PUBLIC_CHARM_D0_PT5_Y0
    print(f"local={local:.12e}")
    print(f"public_D0={PUBLIC_CHARM_D0_PT5_Y0:.12e}")
    print(f"rel_diff={rel:+.6%}   (direct BCFY pseudoscalar only; the feeddown calibration absorbs this)")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("check", choices=("gate", "public-points", "charm-point"))
    gen.add_workspace_options(ap)
    ap.add_argument("--quark", action="append", choices=sorted(gen.QUARKS),
                    help="gate: quark to check; may be repeated (default: bottom and charm)")
    ap.add_argument("--reference-dir", type=Path, default=REFERENCE_DIR,
                    help="gate: directory holding the tracked central grids")
    ap.add_argument("--rtol", type=float, default=1e-3, help="gate: max relative difference allowed")
    ap.add_argument("--grid-workers", type=int, default=1, help="gate: parallel rapidity chunks")
    ap.add_argument("--tolerance", type=float, default=0.05, help="public-points: max relative difference allowed")
    args = ap.parse_args(argv)

    if args.check == "gate":
        out = args.out_dir or Path(tempfile.mkdtemp(prefix="fonll_gate_"))
        ws = gen.Workspace.resolve(args.fonll, out, args.run_dir, args.log_dir)
        ws.out.mkdir(parents=True, exist_ok=True)
        print(f"gate output dir: {ws.out}")
        ok = True
        for quark in args.quark or ["bottom", "charm"]:
            ok = gate_quark(ws, quark, args.reference_dir, args.rtol, args.grid_workers) and ok
        print("GATE PASSED -- central reproduction confirmed" if ok
              else "GATE FAILED -- do not launch the campaign until resolved")
        return 0 if ok else 1
    ws = gen.workspace_from(args)
    if args.check == "public-points":
        return 0 if public_points(ws, args.tolerance) else 1
    charm_point(ws)
    return 0


if __name__ == "__main__":
    sys.exit(main())
