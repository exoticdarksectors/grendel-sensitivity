"""Sequential signal cutflow at one mass point and one coupling.

Rebuilds the acceptance Monte Carlo exactly as the scan does (same seeds,
same chunks, same templates), keeps the full-length selection arrays, weights
every decay by the production weight times the decay density at the given
coupling, and applies the cuts one after another
(``grendel.reco.acceptance.build_cutflow_mc``). The table is therefore the
cutflow of the published scan itself at that point, not a separate sample.

    python -m grendel.cutflow --model hnl --flavor Umu --mass 1.5 --coupling 2e-8
    python -m grendel.cutflow --model bc4 --mass 0.975 --from-results DIR/bc4/sensitivity.csv
    python -m grendel.cutflow --model bc10 --mass 1.23 --coupling 1e-6

``--coupling`` is the scan variable of the model (|U|^2, sin^2 theta, or 1/f
in GeV^-1 for BC10); ``--from-results`` takes the lower island edge of that
mass from a finished scan instead.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

from .geometry.raycast import get_mesh
from .reco.acceptance import build_cutflow_mc
from .scan import MassPoint, ScanConfig, run_point


def cutflow_rows(spec, pt: MassPoint, cfg: ScanConfig, u2: float, mesh=None):
    """The cut-by-cut efficiencies at scan variable ``u2``."""
    mesh = mesh if mesh is not None else get_mesh()
    result = run_point(spec, pt, cfg, mesh, return_mc=True)
    if result is None or result.arrays is None or result.arrays.mc is None:
        raise SystemExit(f"no acceptance Monte Carlo for {pt.tag} (missing inputs?)")
    a = result.arrays
    # scan_u2 integrand at fixed coupling: production weight x path/n x decay density
    lam = a.beta_gamma * a.ctau_ref / u2
    w = a.weights * a.path / a.n_samples
    weights = ((w / lam)[:, None] * np.exp(-a.d / lam[:, None])).reshape(-1)
    if a.sample_w is not None:
        weights = weights * a.sample_w
    return build_cutflow_mc(a.mc, weights), a


def _lower_edge_from_results(path: Path, pt: MassPoint, column: str) -> float:
    with open(path) as fh:
        for row in csv.DictReader(fh):
            if pt.flavor and row.get("flavor") != pt.flavor:
                continue
            if abs(float(row["mass_GeV"]) - pt.mass) < 1e-9:
                return float(row[column])
    raise SystemExit(f"no row for {pt.tag} in {path}")


def main(argv=None) -> int:
    from .models.hnl.scan import resolve_paths
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", choices=("hnl", "bc4", "bc10"), required=True)
    ap.add_argument("--mass", type=float, required=True)
    ap.add_argument("--flavor", default="Umu", help="HNL only")
    ap.add_argument("--coupling", type=float, default=None,
                    help="|U|^2, sin^2 theta, or 1/f [GeV^-1] for the decay-density weights")
    ap.add_argument("--from-results", type=Path, default=None,
                    help="take the lower island edge at this mass from a finished sensitivity.csv")
    ap.add_argument("--decay-samples", type=int, default=None,
                    help="default: the model's scan default")
    ap.add_argument("--templates-dir", type=Path, default=None,
                    help="BC4: cached templates instead of the analytic decay engine")
    ap.add_argument("--vectors-dir", type=Path, default=None)
    ap.add_argument("--geometry-dir", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None, help="analysis directory (paths policy)")
    ap.add_argument("--table", type=Path, default=None, help="write the cutflow CSV here")
    args = ap.parse_args(argv)
    if (args.coupling is None) == (args.from_results is None):
        ap.error("give exactly one of --coupling or --from-results")

    if args.model == "hnl":
        from .models.hnl.spec import HNLSpec
        paths = resolve_paths(args, "hnl")
        spec = HNLSpec(paths)
        pt = MassPoint(args.mass, args.flavor)
        column, to_u2, label = "u2_min", (lambda c: c), "|U|^2"
    elif args.model == "bc4":
        from .models.scalar.spec import ScalarSpec
        paths = resolve_paths(args, "bc4")
        spec = ScalarSpec(paths, templates_dir=args.templates_dir)
        pt = MassPoint(args.mass)
        column, to_u2, label = "u2_min", (lambda c: c), "sin^2 theta"
    else:
        from .models.alp import model as alp_model
        from .models.alp.spec import ALPSpec
        paths = resolve_paths(args, "bc10")
        spec = ALPSpec(paths)
        pt = MassPoint(args.mass)
        # the scan variable is u2 = (1/f / 1/f_ref)^2
        column, to_u2, label = "invf_min", (lambda c: (c / alp_model.INV_F_REF) ** 2), "1/f [GeV^-1]"

    coupling = args.coupling if args.coupling is not None else \
        _lower_edge_from_results(args.from_results, pt, column)
    u2 = to_u2(coupling)
    cfg = ScanConfig(decay_samples=args.decay_samples or spec.decay_samples_default)
    rows, a = cutflow_rows(spec, pt, cfg, u2)

    print(f"{spec.name}  {pt.tag}  {label} = {coupling:.4g}  "
          f"(ctau = {a.ctau_ref / u2:.3g} m)  n_hits = {a.n_hits_eval}")
    print(f"  {'cut':<46}{'cum.eff':>9}{'marg.eff':>10}")
    print("  " + "-" * 65)
    for r in rows:
        print(f"  {r['cut']:<46}{r['efficiency']:>9.4f}{r['marginal_efficiency']:>10.4f}")

    table = args.table or paths.analysis / f"cutflow_{pt.tag.replace('/', '_')}.csv"
    table.parent.mkdir(parents=True, exist_ok=True)
    with open(table, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["cut", "efficiency", "marginal_efficiency"])
        for r in rows:
            wr.writerow([r["cut"], f"{r['efficiency']:.6f}", f"{r['marginal_efficiency']:.6f}"])
    print(f"\nwritten: {table}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
