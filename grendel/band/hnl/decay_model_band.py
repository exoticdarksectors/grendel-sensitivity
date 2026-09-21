"""Decay-model uncertainty band on the HNL exclusion contour.

The decay-model band is ONE coherent nuisance -- the seam-derived
hadronic-width uncertainty ``delta(m)`` (``width_band``) -- driven through
both:

* **lifetime leg:** ``ctau = hbar/Gamma_tot`` -> ``P_decay``; re-scan the
  central run's already-built acceptance Monte Carlo with ``ctau`` rescaled
  by ``1/(1 +/- delta)``. Moves the upper edge ``u2_max``.
* **composition leg:** ``Gamma_had`` shifts the visible/invisible mix ->
  ``vis_frac``; reweight the decay templates by mode from the SAME varied
  widths. Added on top of the lifetime leg; the two share ``Gamma_tot`` and
  partly self-cancel, so they are propagated coherently -- never summed in
  quadrature as independent axes.

Production (four-vectors) and geometry are reused unchanged; reconstruction,
detector response, background and the ``N >= 3`` criterion are held fixed.
The band therefore covers production + decay physics only.

    python -m grendel.band.hnl.decay_model_band --flavor Ue Umu Utau --out band.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from ...geometry.raycast import get_mesh
from ...io.paths import ModelPaths
from ...models.hnl.spec import HNLSpec
from ...scan import MassPoint, ScanConfig, run_point
from . import width_band

# Densified at 3.0-4.4 GeV where the dome closes (the band changes fast there).
# Every entry must be on the HNL mass grid, or the central run has no
# four-vectors for it and the point is dropped (see the warning in run()).
DEFAULT_MASSES = [0.305, 0.5, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8,
                  3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4]


def _f(v):
    return float(v) if v is not None else float("nan")


def run(paths: ModelPaths, flavors, masses, cfg: ScanConfig) -> pd.DataFrame:
    mesh = get_mesh()
    rows = []
    for flavor in flavors:
        for mass in masses:
            d, hf = width_band.delta_and_hadfrac(flavor, mass)
            spec = HNLSpec(paths, width_delta=d, had_frac=hf)
            result = run_point(spec, MassPoint(float(mass), flavor), cfg, mesh)
            if result is None:
                print(f"  {flavor:5s} m={mass}: SKIPPED -- no four-vectors "
                      "(mass not on the grid / not produced in the central run)", flush=True)
                continue
            r = result.row
            if not r.get("has_sensitivity", False):
                continue
            rows.append(r)
            print(f"  {flavor:5s} m={mass:<4} delta={d:.3f} had_frac={hf:.2f}  "
                  f"u2_max={_f(r.get('u2_max')):.3e}  "
                  f"decay-band=[{_f(r.get('u2_max_dm_lo')):.3e}, "
                  f"{_f(r.get('u2_max_dm_hi')):.3e}]", flush=True)
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flavor", nargs="+", default=["Ue", "Umu", "Utau"])
    ap.add_argument("--mass", nargs="+", type=float, default=None)
    ap.add_argument("--decay-samples", type=int, default=50)
    ap.add_argument("--max-hit-events", type=int, default=None,
                    help="optional weighted-resampling cap; omit for exact hits")
    ap.add_argument("--event-chunk", type=int, default=1000,
                    help="exact-hit memory chunk (default: 1000)")
    ap.add_argument("--seed-salt", default="")
    ap.add_argument("--thresholds", nargs="+", type=float, default=[3.0],
                    help="signal-yield thresholds; the first is the exclusion criterion")
    ap.add_argument("--vectors-dir", type=Path, default=None)
    ap.add_argument("--templates-dir", type=Path, default=None)
    ap.add_argument("--analysis-dir", type=Path, default=None,
                    help="the central run's analysis directory (ray-cast cache)")
    ap.add_argument("--geometry-dir", type=Path, default=None)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    if args.decay_samples < 1:
        ap.error("--decay-samples must be >= 1")
    if args.max_hit_events is not None and args.max_hit_events < 1:
        ap.error("--max-hit-events must be >= 1 when provided")
    if args.event_chunk is not None and args.event_chunk < 1:
        ap.error("--event-chunk must be >= 1 when provided")

    paths = ModelPaths.resolve("hnl", vectors=args.vectors_dir, templates=args.templates_dir,
                               analysis=args.analysis_dir, geometry=args.geometry_dir)
    cfg = ScanConfig(decay_samples=args.decay_samples, thresholds=tuple(args.thresholds),
                     max_hit_events=args.max_hit_events, event_chunk=args.event_chunk,
                     seed_salt=args.seed_salt)
    df = run(paths, args.flavor, args.mass or DEFAULT_MASSES, cfg)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"wrote {args.out} ({len(df)} rows with sensitivity)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
