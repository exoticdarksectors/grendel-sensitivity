"""Charged-kaon transport band on the HNL lower edge.

The kaon channel's escape distance (``KAON_D_ESC`` = 1.5 m, varied over
``KAON_D_ESC_RANGE`` = [1, 3] m) is a variation of ``generate_kaon_csvs``:
regenerate the ``Kmeson`` channel with ``--d-esc``, recombine and rescan
once per distance. This assembles the three ``sensitivity.csv`` files into
``kaon_desc_band.csv``, the table ``grendel.models.hnl.plot_money`` draws:
``flavor``, ``mass_GeV`` and each run's lower edge ``u2_min_desc1p0`` /
``u2_min_desc1p5`` / ``u2_min_desc3p0`` (NaN where that run has no
sensitivity or its lower edge is open).

The 1 m and 3 m runs may cover only the kaon-relevant masses, but must cover
the same points, all present in the nominal run.

    python -m grendel.band.hnl.kaon_band --desc1p0 D1/sensitivity.csv \\
        --desc1p5 D15/sensitivity.csv --desc3p0 D3/sensitivity.csv \\
        --out RESULTS/hnl/band/kaon_desc_band.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from ...production.constants import KAON_D_ESC, KAON_D_ESC_RANGE

DISTANCES = (KAON_D_ESC_RANGE[0], KAON_D_ESC, KAON_D_ESC_RANGE[1])   # m: 1.0, 1.5, 3.0


def tag(d_esc: float) -> str:
    return f"desc{d_esc:.1f}".replace(".", "p")


def lower_edge(path, d_esc: float) -> pd.Series:
    """One run's u2_min indexed by (flavor, mass_GeV); NaN unless a closed edge."""
    df = pd.read_csv(path)
    closed = df["has_sensitivity"].fillna(False).astype(bool)
    if "u2_min_open" in df:
        closed &= ~df["u2_min_open"].fillna(False).astype(bool)
    edge = df["u2_min"].where(closed)
    edge.index = pd.MultiIndex.from_frame(df[["flavor", "mass_GeV"]])
    if edge.index.has_duplicates:
        raise ValueError(f"{path}: repeated (flavor, mass_GeV) rows")
    return edge.rename(f"u2_min_{tag(d_esc)}")


def assemble(runs: dict) -> pd.DataFrame:
    """The band table from ``{d_esc: sensitivity.csv}`` of the three runs."""
    edges = {d: lower_edge(runs[d], d) for d in DISTANCES}
    points = edges[DISTANCES[0]].index.sort_values()
    if not points.equals(edges[DISTANCES[-1]].index.sort_values()):
        raise ValueError("the d_esc = 1 m and 3 m runs cover different (flavor, mass) points")
    missing = points.difference(edges[KAON_D_ESC].index)
    if len(missing):
        raise ValueError(f"the nominal d_esc run lacks {len(missing)} of the varied points, "
                         f"e.g. {missing[0]}")
    band = pd.concat([edges[d].reindex(points) for d in DISTANCES], axis=1)
    return band.reset_index()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    for d in DISTANCES:
        ap.add_argument(f"--{tag(d)}", type=Path, required=True,
                        help=f"sensitivity.csv of the d_esc = {d:g} m run")
    ap.add_argument("--out", type=Path, required=True,
                    help="band CSV (plot_money reads <results>/hnl/band/kaon_desc_band.csv)")
    args = ap.parse_args(argv)
    df = assemble({d: getattr(args, tag(d)) for d in DISTANCES})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"wrote {args.out} ({len(df)} rows, d_esc = {', '.join(f'{d:g}' for d in DISTANCES)} m)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
