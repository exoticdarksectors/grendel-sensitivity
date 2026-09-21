"""Draw the figures.

    python -m curves.plot paper --grendel-dir RESULTS --out FIGURES [--benchmark bc4 bc10 hnl]
    python -m curves.plot talk  --grendel-dir RESULTS --out FIGURES [--benchmark ...] [--with-envelope]
    python -m curves.plot diagnostic --grendel-dir RESULTS --out FIGURES

``--grendel-dir`` (or ``GRENDEL_RESULTS_DIR``) holds ``<model>/sensitivity.csv``
and, for the second threshold, ``<model>/sensitivity_nsig10.csv``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .results import MODELS, GrendelResults


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m curves.plot", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("layout", choices=("paper", "talk", "diagnostic"))
    ap.add_argument("--grendel-dir", type=Path, default=None, help="results tree (default: $GRENDEL_RESULTS_DIR)")
    ap.add_argument("--out", type=Path, default=Path("figures"), help="figure directory (default: ./figures)")
    ap.add_argument("--benchmark", nargs="+", choices=MODELS, default=list(MODELS))
    ap.add_argument("--secondary-threshold", type=float, default=10.0,
                    help="second signal threshold drawn dashed inside the N_sig >= 3 band; 0 to omit")
    ap.add_argument("--suffix", default="", help="appended to the BC4/BC10 figure stems")
    ap.add_argument("--stem", default="hnlimits_grendel_paper", help="the HNL figure stem")
    ap.add_argument("--with-envelope", action="store_true",
                    help="talk layout: overlay the single-source variation envelope (needs <model>/band/)")
    args = ap.parse_args(argv)

    results = GrendelResults(args.grendel_dir)
    secondary = args.secondary_threshold or None
    if args.layout == "paper":
        from .paper import render_bc_paper, render_hnl_paper
        for benchmark in args.benchmark:
            if benchmark == "hnl":
                print(render_hnl_paper(results, args.out, stem=args.stem, secondary_threshold=secondary))
            else:
                print(render_bc_paper(benchmark, results, args.out, secondary_threshold=secondary,
                                      output_suffix=args.suffix))
    elif args.layout == "talk":
        from .talk import render_bc_talk, render_hnl_talk
        for benchmark in args.benchmark:
            if benchmark == "hnl":
                for out in render_hnl_talk(results, args.out, secondary_threshold=secondary):
                    print(out)
            else:
                print(render_bc_talk(benchmark, results, args.out, with_envelope=args.with_envelope))
    else:
        from .diagnostic import render_diagnostics
        for out in render_diagnostics(results, args.out, args.benchmark):
            print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
