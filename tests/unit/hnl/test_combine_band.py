"""Topology-accounting tests for FONLL band combination."""

import json

import numpy as np
import pandas as pd

from grendel.band.quadrature import combine_band


def _write_curve(path, *, has_sensitivity=True, u2_min=1e-8,
                 u2_max=1e-4, u2_max_open=False):
    pd.DataFrame([{
        "flavor": "Ue",
        "mass_GeV": 1.0,
        "has_sensitivity": has_sensitivity,
        "u2_min": u2_min,
        "u2_max": u2_max,
        "u2_min_open": False,
        "u2_max_open": u2_max_open,
    }]).to_csv(path, index=False)


def _registry(tmp_path, varied_kwargs):
    central = tmp_path / "central.csv"
    varied = tmp_path / "varied.csv"
    _write_curve(central)
    _write_curve(varied, **varied_kwargs)
    registry = tmp_path / "registry.json"
    registry.write_text(json.dumps({
        "variations": [
            {"name": "central", "axis": "central",
             "sensitivity_csv": str(central)},
            {"name": "scale_test", "axis": "scale",
             "sensitivity_csv": str(varied)},
        ]
    }))
    return registry


def test_combine_band_records_stable_member_counts(tmp_path):
    row = combine_band(_registry(tmp_path, {"u2_min": 2e-8})).iloc[0]

    assert not bool(row.u2_min_topology_changed)
    assert row.u2_min_n_members_expected == 1
    assert row.u2_min_n_members_sensitive == 1
    assert row.u2_min_n_members_finite == 1
    assert np.isfinite(row.u2_min_band_lo)
    assert np.isfinite(row.u2_min_band_hi)


def test_combine_band_suppresses_ribbon_on_lost_island(tmp_path):
    row = combine_band(_registry(tmp_path, {
        "has_sensitivity": False,
        "u2_min": np.nan,
        "u2_max": np.nan,
    })).iloc[0]

    assert bool(row.u2_min_topology_changed)
    assert row.u2_min_n_members_sensitive == 0
    assert np.isnan(row.u2_min_band_lo)
    assert np.isnan(row.u2_min_band_hi)
