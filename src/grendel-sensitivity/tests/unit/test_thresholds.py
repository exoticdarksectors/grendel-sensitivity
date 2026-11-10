import pandas as pd
import pytest

from grendel.io.thresholds import split_threshold


def test_split_moves_the_secondary_island_into_place():
    frame = pd.DataFrame([{
        "u2_min": 1e-9, "u2_max": 1e-6, "u2_min_open": False, "u2_max_open": False, "peak_N": 40.0,
        "peak_u2": 1e-8, "has_sensitivity": True,
        "u2_min_N10": 3e-9, "u2_max_N10": 3e-7, "peak_u2_N10": 1e-8, "has_sensitivity_N10": True,
        "u2_min_open_N10": False, "u2_max_open_N10": False, "mass_GeV": 1.0, "n_events": 5}])
    out = split_threshold(frame, 10.0)
    assert list(out.columns) == ["u2_min", "u2_max", "u2_min_open", "u2_max_open", "peak_N", "peak_u2",
                                 "has_sensitivity", "mass_GeV", "n_events"]
    assert out.loc[0, "u2_min"] == 3e-9 and out.loc[0, "u2_max"] == 3e-7
    assert out.loc[0, "peak_N"] == 40.0      # the yield curve is the same
    with pytest.raises(ValueError, match="not solved"):
        split_threshold(frame, 5.0)
