"""Closing the BC4/BC10 islands where the yield crosses the threshold."""
import numpy as np
import pytest

from curves.plot.bc_panel import insert_threshold_tips, sensitive_segments
from curves.plot.common import SIGNAL_THRESHOLD


def _row(mass, peak, coupling, sensitive):
    return {"mass_GeV": mass, "c_min": coupling / 2 if sensitive else np.nan,
            "c_max": coupling * 2 if sensitive else np.nan, "peak_N": peak, "peak_c": coupling,
            "has_sensitivity": sensitive}


def test_every_supported_sensitivity_transition_gets_a_closed_tip():
    rows = [_row(1.0, 4.0, 1.0e-8, True), _row(2.0, 2.0, 2.0e-8, False),
            _row(3.0, 1.5, 3.0e-8, False), _row(4.0, 6.0, 6.0e-8, True)]
    tips = [row for row in insert_threshold_tips(rows) if row["peak_N"] == SIGNAL_THRESHOLD]
    assert len(tips) == 2
    assert 1.0 < tips[0]["mass_GeV"] < 2.0
    assert 3.0 < tips[1]["mass_GeV"] < 4.0
    for tip in tips:
        assert tip["has_sensitivity"]
        assert tip["c_min"] == pytest.approx(tip["peak_c"])
        assert tip["c_max"] == pytest.approx(tip["peak_c"])


def test_unsupported_nonfinite_row_is_not_bridged():
    rows = [_row(0.52, 20.0, 1.0e-7, True), _row(0.54, np.nan, np.nan, False), _row(0.56, 18.0, 1.2e-7, True)]
    assert insert_threshold_tips(rows) == rows


def test_final_closure_remains_interpolated():
    rows = [_row(3.2, 3.8, 2.2e-8, True), _row(3.3, 2.9, 2.4e-8, False)]
    result = insert_threshold_tips(rows)
    assert len(result) == 3
    assert result[1]["peak_N"] == SIGNAL_THRESHOLD
    assert 3.2 < result[1]["mass_GeV"] < 3.3


def test_tips_close_at_the_requested_threshold():
    # Peak yield falls 40 -> 4 across the step, so the N = 10 crossing sits
    # strictly inside it while the N = 3 crossing never occurs.
    rows = [_row(1.0, 40.0, 1.0e-8, True), _row(2.0, 4.0, 2.0e-8, False)]
    tips = [row for row in insert_threshold_tips(rows, threshold=10.0) if row["peak_N"] == 10.0]
    assert len(tips) == 1
    assert 1.0 < tips[0]["mass_GeV"] < 2.0
    assert not [row for row in insert_threshold_tips(rows) if row["peak_N"] == SIGNAL_THRESHOLD]


def test_default_threshold_is_unchanged():
    rows = [_row(1.0, 4.0, 1.0e-8, True), _row(2.0, 2.0, 2.0e-8, False)]
    assert insert_threshold_tips(rows) == insert_threshold_tips(rows, threshold=SIGNAL_THRESHOLD)


def test_narrow_islands_are_dropped_only_when_asked(capsys):
    rows = [_row(1.0, 5.0, 1e-8, True), _row(1.1, 5.0, 1e-8, True), _row(1.2, 2.0, 1e-8, False),
            _row(1.40, 3.19, 1e-8, True), _row(1.41, 3.16, 1e-8, True), _row(1.5, 2.0, 1e-8, False)]
    assert [len(s) for s in sensitive_segments(rows)] == [2, 2]
    assert [len(s) for s in sensitive_segments(rows, min_width_GeV=0.05)] == [2]
    assert "dropped island 1.4000-1.4100 GeV" in capsys.readouterr().err
