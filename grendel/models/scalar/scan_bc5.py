"""Scan the BC5 dark-scalar sensitivity: the (m_S, sin^2 theta) island at BR(h -> SS) = 0.01 over the
whole mass range, to m_h / 2."""
from __future__ import annotations

import argparse
import sys

from ...scan import NoResultsError, ScanConfig, run_scan
from ..hnl.scan import add_common_options, resolve_paths
from . import production_bc5 as production
from .spec_bc5 import QuarticScalarSpec


def produce_missing(masses, vectors_dir, n_pool, n_higgs, seed, *, br_hss, quartic_b,
                    force=False, high_pt_tilt_scale=None, nominal_mixture_fraction=0.5) -> int:
    """Write the CSVs that are missing (all with ``force``) from one shared b-hadron pool and one
    shared Higgs sample."""
    vectors_dir.mkdir(parents=True, exist_ok=True)
    if force:
        to_produce = list(masses)
    else:
        to_produce = [m for m in masses
                      if not (vectors_dir / f"mS_{production.format_mass_for_filename(m)}.csv").exists()]
    if not to_produce:
        print(f"reusing existing four-vector CSVs for all {len(masses)} masses "
              "(--force-produce to regenerate)")
        return 0
    rngs = production.make_streams(seed)
    pool, higgs, sigma_bottom = production.shared_pools(
        rngs, n_pool, n_higgs, br_hss=br_hss, high_pt_tilt_scale=high_pt_tilt_scale,
        nominal_mixture_fraction=nominal_mixture_fraction)
    sigma_h = higgs["sigma_pb"] if higgs is not None else 0.0
    print(f"sigma_FONLL(bottom) = {sigma_bottom:.3e} pb; sigma(pp -> h) = {sigma_h:.2f} pb x "
          f"BR(h -> SS) = {br_hss:g}; quartic B decays {'on' if quartic_b else 'off'}; "
          f"producing {len(to_produce)}/{len(masses)} masses")
    for m_S in to_produce:
        path, counts = production.write_bc5_csv(
            m_S, vectors_dir, rngs, pool=pool, higgs=higgs, sigma_bottom=sigma_bottom,
            br_hss=br_hss, quartic_b=quartic_b)
        print(f"  m_S={m_S:.3f}: mixing {counts['mixing']}, h->SS {counts['hSS']}, "
              f"B->X SS {counts['BSS']} rows -> {path.name}")
    return len(to_produce)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add_common_options(parser, "bc5", QuarticScalarSpec.decay_samples_default)
    parser.add_argument("--n-pool", type=int, default=production.production.N_POOL_DEFAULT,
                        help="parent b-hadron pool size when producing four-vectors")
    parser.add_argument("--n-higgs", type=int, default=production.N_HIGGS_DEFAULT,
                        help="Higgs sample size when producing four-vectors")
    parser.add_argument("--seed", type=int, default=42, help="production seed")
    parser.add_argument("--br-hss", type=float, default=production.model.BR_HSS_BC5,
                        help="BR(h -> SS) fixing the quartic coupling (0: no quartic channels)")
    parser.add_argument("--no-quartic-b", action="store_true",
                        help="drop b -> X_s S S and B_s -> S S (cross-checks)")
    parser.add_argument("--force-produce", action="store_true",
                        help="regenerate the four-vector CSVs even when present")
    parser.add_argument("--high-pt-tilt-scale", type=float, default=None,
                        help="importance-sample a nominal/high-pT FONLL mixture tilted by exp(pT/scale)")
    parser.add_argument("--nominal-mixture-fraction", type=float, default=0.5,
                        help="nominal fraction of the optional high-pT proposal mixture")
    parser.add_argument("--width-scheme", default="winkler",
                        help="width model of the analytic decay engine")
    parser.add_argument("--seed-offset", type=int, default=0,
                        help="added to every per-mass reconstruction seed")
    args = parser.parse_args(argv)
    if args.decay_samples < 1:
        parser.error("--decay-samples must be >= 1")

    paths = resolve_paths(args, "bc5")
    spec = QuarticScalarSpec(paths, templates_dir=args.templates_dir,
                             width_scheme=args.width_scheme, br_hss=args.br_hss)
    masses = args.mass or production.MASS_GRID_BC5
    produce_missing(masses, paths.vectors, args.n_pool, args.n_higgs, args.seed,
                    br_hss=args.br_hss, quartic_b=not args.no_quartic_b,
                    force=args.force_produce, high_pt_tilt_scale=args.high_pt_tilt_scale,
                    nominal_mixture_fraction=args.nominal_mixture_fraction)
    cfg = ScanConfig(decay_samples=args.decay_samples, thresholds=tuple(args.thresholds),
                     seed_offset=args.seed_offset, force_geometry=args.force_geometry)
    print(paths.describe())
    try:
        run_scan(spec, spec.points(masses), cfg, paths.analysis,
                 workers=args.workers, resume=args.resume)
    except NoResultsError as err:
        print(f"error: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
