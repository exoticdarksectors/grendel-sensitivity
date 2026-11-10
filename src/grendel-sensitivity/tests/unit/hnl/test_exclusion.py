"""Tests for resolved and scan-limited exclusion boundaries."""

import numpy as np
import pandas as pd

from grendel.reco.exclusion import find_exclusion_band, find_exclusion_band_refined
from grendel.models.hnl.plot_exclusion import _sensitive_segments, _thin_marker_mask


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


def test_refined_band_calls_exact_evaluator_at_edges_and_peak():
    # N(x) = A x^2 exp(-B x) peaks at x=2/B.  A deliberately coarse scan
    # makes interpolation visibly inadequate while still bracketing the island.
    amplitude = 40.0
    slope = 3.0

    def signal(x):
        return amplitude * x**2 * np.exp(-slope * x)

    u2 = np.logspace(-2.0, 1.0, 18)
    result = find_exclusion_band_refined(
        u2, signal(u2), signal, N_threshold=1.0,
    )

    assert result["has_sensitivity"]
    assert np.isclose(result["peak_u2"], 2.0 / slope, rtol=2.0e-5)
    assert np.isclose(result["peak_N"], signal(2.0 / slope), rtol=1.0e-9)
    assert np.isclose(signal(result["u2_min"]), 1.0, rtol=2.0e-5)
    assert np.isclose(signal(result["u2_max"]), 1.0, rtol=2.0e-5)


def test_open_edge_markers_thin_to_sparse_evenly_spaced_subset():
    open_mask = np.ones(40, dtype=bool)

    thinned = _thin_marker_mask(open_mask, max_markers=12)

    assert thinned.sum() == 12
    assert thinned[0] and thinned[-1]          # keep the span endpoints
    # a short fully-open span is left untouched (nothing to thin)
    short = np.ones(5, dtype=bool)
    assert _thin_marker_mask(short, max_markers=12) is short


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
