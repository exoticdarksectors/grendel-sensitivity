"""Scan the BC4 dark-scalar sensitivity: the (m_S, sin^2 theta) island.

    python -m grendel.models.scalar.scan                   # full grid; produces only
                                                           # missing four-vectors
    python -m grendel.models.scalar.scan --force-produce   # regenerate the four-vectors
    python -m grendel.models.scalar.scan --mass 0.5 1 2 --n-pool 100000
    python -m grendel.models.scalar.scan --templates-dir DIR   # exHad decays
    python -m grendel.models.scalar.scan --thresholds 3 10

Four-vectors are generated for masses that have none (or for all with
``--force-produce``); untouched files keep their mtime so the ray-cast cache
stays valid and a re-scan only redoes the decay Monte Carlo, reconstruction
and coupling scan. Paths follow ``grendel.io.paths`` (``GRENDEL_BC4_*``).
"""
from __future__ import annotations

import argparse
import sys

import numpy as np

from ...production.fonll.fonll_parser import get_sigma_total
from ...production.fonll.meson_sampler import sample_meson_4vectors
from ...scan import ScanConfig, run_scan
from ..hnl.scan import add_common_options, resolve_paths
from . import production
from .spec import ScalarSpec


def produce_missing(masses, vectors_dir, n_pool, seed, *, force=False,
                    high_pt_tilt_scale=None, nominal_mixture_fraction=0.5) -> int:
    """Write the four-vector CSVs that are missing (all of them with
    ``force``), from one shared parent pool. Returns how many were written."""
    vectors_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    if force:
        to_produce = list(masses)
    else:
        to_produce = [m for m in masses
                      if not (vectors_dir / f"mS_{production._mass_label(m)}.csv").exists()]
    if not to_produce:
        print(f"reusing existing four-vector CSVs for all {len(masses)} masses "
              "(--force-produce to regenerate)")
        return 0
    sigma_bottom = get_sigma_total("bottom")
    pool = sample_meson_4vectors(n_pool, "bottom", rng=rng,
                                 high_pt_tilt_scale=high_pt_tilt_scale,
                                 nominal_mixture_fraction=nominal_mixture_fraction)
    print(f"sigma_FONLL(bottom) = {sigma_bottom:.3e} pb; "
          f"producing {len(to_produce)}/{len(masses)} masses")
    for m_S in to_produce:
        production.write_scalar_csv(m_S, vectors_dir, n_pool, rng,
                                    sigma_bottom=sigma_bottom, pool=pool)
    return len(to_produce)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add_common_options(parser, "bc4", ScalarSpec.decay_samples_default)
    parser.add_argument("--n-pool", type=int, default=production.N_POOL_DEFAULT,
                        help="parent b-hadron pool size when producing four-vectors")
    parser.add_argument("--seed", type=int, default=42, help="production seed")
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

    paths = resolve_paths(args, "bc4")
    spec = ScalarSpec(paths, templates_dir=args.templates_dir, width_scheme=args.width_scheme)
    out_csv = paths.analysis / "sensitivity.csv"
    if args.plot_only:
        if not out_csv.exists():
            print(f"no results to plot: {out_csv}")
            return 1
        spec.plot(out_csv, paths.analysis)
        return 0

    masses = args.mass or production.MASS_GRID
    produce_missing(masses, paths.vectors, args.n_pool, args.seed,
                    force=args.force_produce,
                    high_pt_tilt_scale=args.high_pt_tilt_scale,
                    nominal_mixture_fraction=args.nominal_mixture_fraction)
    cfg = ScanConfig(decay_samples=args.decay_samples, thresholds=tuple(args.thresholds),
                     seed_offset=args.seed_offset, force_geometry=args.force_geometry)
    print(paths.describe())
    run_scan(spec, spec.points(masses), cfg, paths.analysis,
             workers=args.workers, resume=args.resume)
    return 0


if __name__ == "__main__":
    sys.exit(main())
