"""Tests for resolved and scan-limited exclusion boundaries."""

import numpy as np
import pandas as pd

from analysis.exclusion import find_exclusion_band
from analysis.plot_exclusion import _sensitive_segments


def test_closed_exclusion_band_has_two_resolved_edges():
    u2 = np.logspace(-5, -2, 4)
    n_signal = np.array([1.0, 4.0, 5.0, 1.0])

    result = find_exclusion_band(u2, n_signal, N_threshold=3.0)

    assert result["has_sensitivity"]
    assert not result["u2_min_open"]
    assert not result["u2_max_open"]
    assert np.isfinite(result["u2_min"])
    assert np.isfinite(result["u2_max"])


def test_upper_edge_is_open_when_scan_ends_above_threshold():
    u2 = np.logspace(-5, -2, 4)
    n_signal = np.array([1.0, 4.0, 5.0, 6.0])

    result = find_exclusion_band(u2, n_signal, N_threshold=3.0)

    assert result["has_sensitivity"]
    assert not result["u2_min_open"]
    assert result["u2_max_open"]
    assert np.isfinite(result["u2_min"])
    assert np.isnan(result["u2_max"])


def test_lower_edge_is_open_when_scan_starts_above_threshold():
    u2 = np.logspace(-5, -2, 4)
    n_signal = np.array([6.0, 5.0, 4.0, 1.0])

    result = find_exclusion_band(u2, n_signal, N_threshold=3.0)

    assert result["has_sensitivity"]
    assert result["u2_min_open"]
    assert not result["u2_max_open"]
    assert np.isnan(result["u2_min"])
    assert np.isfinite(result["u2_max"])


def test_no_sensitivity_has_no_open_edges():
    u2 = np.logspace(-5, -2, 4)
    n_signal = np.array([0.0, 1.0, 2.0, 1.0])

    result = find_exclusion_band(u2, n_signal, N_threshold=3.0)

    assert not result["has_sensitivity"]
    assert not result["u2_min_open"]
    assert not result["u2_max_open"]
    assert np.isnan(result["u2_min"])
    assert np.isnan(result["u2_max"])


def test_plot_segments_do_not_bridge_insensitive_mass_points():
    rows = pd.DataFrame({
        "mass_GeV": [0.2, 0.3, 0.4, 0.5],
        "has_sensitivity": [True, True, False, True],
    })

    segments = list(_sensitive_segments(rows))

    assert [segment["mass_GeV"].tolist() for segment in segments] == [
        [0.2, 0.3],
        [0.5],
    ]
