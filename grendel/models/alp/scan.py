"""Scan the BC10 fermiophilic-ALP sensitivity: the (m_a, 1/f) island.

    python -m grendel.models.alp.scan                     # full grid
    python -m grendel.models.alp.scan --mass 1.0 2.0      # a mass subset
    python -m grendel.models.alp.scan --thresholds 3 10   # also solve N >= 10

Four-vectors come from ``python -m grendel.models.alp.production`` and the
decay templates from the template generators; paths follow
``grendel.io.paths`` (``GRENDEL_BC10_*``). Results land in
``<out>/sensitivity.csv``.
"""
from __future__ import annotations

import argparse
import sys

from ...scan import ScanConfig, run_scan
from ..hnl.scan import add_common_options, resolve_paths
from .mass_grid import ALP_MASS_GRID
from .spec import ALPSpec


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add_common_options(parser, "bc10", ALPSpec.decay_samples_default)
    parser.add_argument("--reco-seed-offset", type=int, default=0,
                        help="added to every deterministic per-mass reconstruction seed")
    args = parser.parse_args(argv)
    if args.decay_samples <= 0:
        parser.error("--decay-samples must be positive")

    paths = resolve_paths(args, "bc10")
    spec = ALPSpec(paths)
    out_csv = paths.analysis / "sensitivity.csv"
    if args.plot_only:
        if not out_csv.exists():
            print(f"no results to plot: {out_csv}")
            return 1
        spec.plot(out_csv, paths.analysis)
        return 0

    masses = args.mass or ALP_MASS_GRID
    cfg = ScanConfig(decay_samples=args.decay_samples, thresholds=tuple(args.thresholds),
                     seed_offset=args.reco_seed_offset, force_geometry=args.force_geometry)
    print(paths.describe())
    run_scan(spec, spec.points(masses), cfg, paths.analysis,
             workers=args.workers, resume=args.resume)
    return 0


if __name__ == "__main__":
    sys.exit(main())
