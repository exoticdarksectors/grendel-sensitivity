"""Rebuild the HNL band's member curves from the flat raw table.

A campaign leaves one sensitivity CSV per variation and a registry that
points at them. The compact record of a campaign is instead one table with
every variation's rows stacked (``variation_name``, ``variation_axis``,
``run_tag``, ... plus the sensitivity columns). This module turns that table
back into the per-variation curves :func:`grendel.band.quadrature.combine_curves`
consumes, so a band can be recombined from the raw table alone.
"""
from __future__ import annotations

import pandas as pd

from ..quadrature import AXES, combine_curves

META_COLUMNS = ("variation_index", "variation_name", "variation_axis",
                "variation_direction", "run_tag", "bottom_grid_tag", "charm_grid_tag")


def curves_from_raw(raw: pd.DataFrame) -> tuple[dict[str, pd.DataFrame], dict[str, list[str]]]:
    """(curves by variation name indexed by (flavor, mass), axis -> names)."""
    curves = {}
    axes: dict[str, list[str]] = {axis: [] for axis in AXES}
    drop = [c for c in META_COLUMNS if c in raw.columns]
    for name, df in raw.groupby("variation_name", sort=False):
        axes[str(df["variation_axis"].iloc[0])].append(str(name))
        curves[str(name)] = df.drop(columns=drop).set_index(["flavor", "mass_GeV"]).sort_index()
    return curves, axes


def combine_raw(raw: pd.DataFrame, *, alphas=None, suppress=frozenset()) -> pd.DataFrame:
    curves, axes = curves_from_raw(raw)
    return combine_curves(curves, axes, alphas=alphas, suppress=suppress)


def raw_from_registry(registry: dict) -> pd.DataFrame:
    """The flat table of a campaign registry (the inverse of :func:`curves_from_raw`)."""
    frames = []
    for index, v in enumerate(registry["variations"]):
        df = pd.read_csv(v["sensitivity_csv"])
        meta = {"variation_index": index, "variation_name": v["name"],
                "variation_axis": v["axis"], "variation_direction": v.get("direction", ""),
                "run_tag": v.get("run_tag", ""), "bottom_grid_tag": v.get("bottom_grid_tag", ""),
                "charm_grid_tag": v.get("charm_grid_tag", "")}
        for column, value in reversed(list(meta.items())):
            df.insert(0, column, value)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)
