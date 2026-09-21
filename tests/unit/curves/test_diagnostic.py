import numpy as np
import pandas as pd

from curves.plot.diagnostic import bc_vertical_edges, break_jumps, dex_densify


def test_dex_densify_does_not_bridge_or_extrapolate_unsupported_anchors():
    masses = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    central = np.ones(5)
    lower = np.array([0.8, 0.8, np.nan, 0.7, 0.7])
    upper = np.array([1.2, 1.2, np.nan, 1.3, 1.3])
    dense = np.arange(0.5, 5.6, 0.5)
    _, dense_lower, dense_upper = dex_densify(masses, central, lower, upper, dense, np.ones(len(dense)))
    supported = np.isin(dense, [1.0, 1.5, 2.0, 4.0, 4.5, 5.0])
    assert np.isfinite(dense_lower[supported]).all()
    assert np.isfinite(dense_upper[supported]).all()
    assert np.isnan(dense_lower[~supported]).all()
    assert np.isnan(dense_upper[~supported]).all()


def test_bc_vertical_edges_follow_plotted_strength_not_suffix_name():
    lower, upper = bc_vertical_edges(pd.DataFrame({"u2_min_bc_lo": [1.2e-8], "u2_min_bc_hi": [0.8e-8]}))
    assert lower[0] == 0.8e-8
    assert upper[0] == 1.2e-8


def test_break_jumps_inserts_a_gap_at_large_steps():
    x, y = break_jumps([1.0, 2.0, 3.0], [1e-8, 1e-3, 2e-3])
    assert len(x) == 4 and np.isnan(x[1]) and np.isnan(y[1])
