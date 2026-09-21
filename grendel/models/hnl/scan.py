"""Scan the HNL sensitivity: (m_N, |U|^2) exclusion band per flavor.

    python -m grendel.models.hnl.scan                         # Umu, full grid
    python -m grendel.models.hnl.scan --flavor Ue Umu Utau    # every flavor
    python -m grendel.models.hnl.scan --mass 0.5 1.0 2.0      # a mass subset
    python -m grendel.models.hnl.scan --thresholds 3 10       # also solve N >= 10
    python -m grendel.models.hnl.scan --workers 6

Inputs and outputs follow ``grendel.io.paths`` (``GRENDEL_HNL_*`` or the
``--vectors-dir``/``--templates-dir``/``--out``/``--geometry-dir`` options);
results land in ``<out>/sensitivity.csv``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from ...constants import FLAVORS
from ...io.paths import ModelPaths
from ...scan import ScanConfig, run_scan
from .mass_grid import ANALYSIS_MASS_MAX, MASS_GRID
from .spec import HNLSpec


def add_common_options(parser: argparse.ArgumentParser, model: str, decay_samples: int) -> None:
    key = model.upper()
    parser.add_argument("--mass", nargs="+", type=float, default=None,
                        help="masses in GeV (default: the model's grid)")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--decay-samples", type=int, default=decay_samples,
                        help=f"decay positions sampled per hitting event (default {decay_samples})")
    parser.add_argument("--thresholds", nargs="+", type=float, default=[3.0],
                        help="signal-yield thresholds; the first is the exclusion "
                             "criterion, the rest are solved on the same yield curve "
                             "(default: 3)")
    parser.add_argument("--force-geometry", action="store_true",
                        help="recompute the ray-cast cache")
    parser.add_argument("--resume", action="store_true",
                        help="keep rows already in <out>/sensitivity.csv and skip their points")
    parser.add_argument("--plot-only", action="store_true",
                        help="re-plot from the existing sensitivity.csv")
    parser.add_argument("--vectors-dir", type=Path, default=None,
                        help=f"four-vector CSVs (default: $GRENDEL_{key}_VECTORS_DIR)")
    parser.add_argument("--templates-dir", type=Path, default=None,
                        help=f"decay templates (default: $GRENDEL_{key}_TEMPLATES_DIR)")
    parser.add_argument("--out", type=Path, default=None,
                        help=f"analysis directory (default: $GRENDEL_{key}_ANALYSIS_DIR)")
    parser.add_argument("--geometry-dir", type=Path, default=None,
                        help=f"ray-cast cache (default: <out>/geometry_cache or $GRENDEL_{key}_GEOMETRY_DIR)")


def resolve_paths(args, model: str) -> ModelPaths:
    return ModelPaths.resolve(model, vectors=args.vectors_dir, templates=args.templates_dir,
                              analysis=args.out, geometry=args.geometry_dir)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--flavor", nargs="+", default=None, choices=FLAVORS,
                        help=f"flavors to scan (default: {list(HNLSpec.default_flavors)})")
    add_common_options(parser, "hnl", HNLSpec.decay_samples_default)
    parser.add_argument("--max-hit-events", type=int, default=None,
                        help="approximate mode: weighted-resample at most this many hitting events per point")
    parser.add_argument("--event-chunk", type=int, default=None,
                        help="bound peak memory by reconstructing this many hitting events at a time")
    parser.add_argument("--seed-salt", default="",
                        help="deterministic salt for independent numerical-control repeats")
    args = parser.parse_args(argv)

    if args.decay_samples < 1:
        parser.error("--decay-samples must be >= 1")
    if args.max_hit_events is not None and args.max_hit_events < 1:
        parser.error("--max-hit-events must be >= 1 when provided")
    if args.event_chunk is not None and args.event_chunk < 1:
        parser.error("--event-chunk must be >= 1 when provided")

    paths = resolve_paths(args, "hnl")
    spec = HNLSpec(paths)
    out_csv = paths.analysis / "sensitivity.csv"
    if args.plot_only:
        if not out_csv.exists():
            print(f"no results to plot: {out_csv}")
            return 1
        spec.plot(out_csv, paths.analysis)
        return 0

    flavors = args.flavor or list(HNLSpec.default_flavors)
    masses = args.mass or [m for m in MASS_GRID if m <= ANALYSIS_MASS_MAX]
    cfg = ScanConfig(decay_samples=args.decay_samples, thresholds=tuple(args.thresholds),
                     max_hit_events=args.max_hit_events, event_chunk=args.event_chunk,
                     seed_salt=args.seed_salt, force_geometry=args.force_geometry)
    print(paths.describe())
    run_scan(spec, spec.points(masses, flavors), cfg, paths.analysis,
             workers=args.workers, resume=args.resume)
    return 0


if __name__ == "__main__":
    sys.exit(main())
