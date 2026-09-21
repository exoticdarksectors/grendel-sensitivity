import numpy as np
import pandas as pd

from grendel.models.scalar.plot_exclusion import SIGNAL_THRESHOLD, _with_threshold_tips


def _row(mass, peak, coupling, sensitive):
    return {
        "mass_GeV": mass,
        "has_sensitivity": sensitive,
        "peak_N": peak,
        "peak_u2": coupling,
        "u2_min": coupling / 2 if sensitive else np.nan,
        "u2_max": coupling * 2 if sensitive else np.nan,
        "u2_min_open": False,
        "u2_max_open": False,
    }


def test_plot_inserts_finite_threshold_tip():
    frame = pd.DataFrame([
        _row(3.7, 5.2, 5.7e-11, True),
        _row(3.8, 2.9, 3.7e-11, False),
    ])
    result = _with_threshold_tips(frame)
    tip = result[result.peak_N == SIGNAL_THRESHOLD].iloc[0]
    assert 3.7 < tip.mass_GeV < 3.8
    assert tip.u2_min == tip.peak_u2 == tip.u2_max


def test_plot_does_not_bridge_nonfinite_transition():
    frame = pd.DataFrame([
        _row(0.52, 10.0, 1.0e-8, True),
        _row(0.54, np.nan, np.nan, False),
        _row(0.56, 9.0, 1.2e-8, True),
    ])
    assert len(_with_threshold_tips(frame)) == len(frame)
