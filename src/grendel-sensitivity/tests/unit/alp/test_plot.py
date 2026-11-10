import numpy as np
import pandas as pd
import pytest

from grendel.models.alp.plot import SIGNAL_THRESHOLD, _insert_threshold_tips


def _row(mass, peak, coupling, sensitive):
    return {
        "mass_GeV": mass,
        "has_sensitivity": sensitive,
        "peak_N": peak,
        "peak_invf": coupling,
        "invf_min": coupling / 2.0 if sensitive else np.nan,
        "invf_max": coupling * 2.0 if sensitive else np.nan,
        "invf_min_open": False,
        "invf_max_open": False,
    }


def test_finite_sensitivity_transitions_get_closed_tips():
    frame = pd.DataFrame(
        [
            _row(1.0, 4.0, 1.0e-8, True),
            _row(2.0, 2.0, 2.0e-8, False),
            _row(3.0, 1.5, 3.0e-8, False),
            _row(4.0, 6.0, 6.0e-8, True),
        ]
    )
    result = _insert_threshold_tips(frame)
    tips = result[result["peak_N"] == SIGNAL_THRESHOLD]

    assert len(tips) == 2
    assert 1.0 < tips.iloc[0]["mass_GeV"] < 2.0
    assert 3.0 < tips.iloc[1]["mass_GeV"] < 4.0
    assert np.allclose(tips["invf_min"], tips["peak_invf"])
    assert np.allclose(tips["invf_max"], tips["peak_invf"])


def test_nonfinite_unsupported_row_is_not_bridged():
    frame = pd.DataFrame(
        [
            _row(0.52, 20.0, 1.0e-7, True),
            _row(0.54, np.nan, np.nan, False),
            _row(0.56, 18.0, 1.2e-7, True),
        ]
    )
    result = _insert_threshold_tips(frame)

    assert len(result) == len(frame)
    pd.testing.assert_frame_equal(
        result.reset_index(drop=True), frame.reset_index(drop=True)
    )


def test_tip_interpolates_yield_and_coupling_logarithmically():
    frame = pd.DataFrame(
        [_row(3.2, 3.8, 2.2e-8, True), _row(3.3, 2.9, 2.4e-8, False)]
    )
    result = _insert_threshold_tips(frame)
    tip = result.iloc[1]

    assert tip["peak_N"] == pytest.approx(SIGNAL_THRESHOLD)
    assert 3.2 < tip["mass_GeV"] < 3.3
    assert tip["invf_min"] == pytest.approx(tip["invf_max"])
