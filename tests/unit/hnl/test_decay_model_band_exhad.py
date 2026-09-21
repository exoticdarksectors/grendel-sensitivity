"""The exHad member of the HNL decay-model band (analysis/decay_model_band_exhad.py)."""
import numpy as np
import pandas as pd
import pytest

from grendel.band.hnl import decay_model_band_exhad as dmx


def _row(flavor, mass, u2_min, u2_max, *, n_hits=100, dm_lo=None, dm_hi=None,
         has_sensitivity=True, u2_max_open=False):
    return {
        "flavor": flavor, "mass_GeV": mass, "has_sensitivity": has_sensitivity,
        "u2_min": u2_min, "u2_max": u2_max, "u2_min_open": False, "u2_max_open": u2_max_open,
        "peak_N": 10.0, "peak_u2": 1e-6, "n_events": 1000, "n_hits": n_hits, "n_hits_eval": n_hits,
        "u2_min_dm_lo": dm_lo[0] if dm_lo else np.nan, "u2_min_dm_hi": dm_hi[0] if dm_hi else np.nan,
        "u2_max_dm_lo": dm_lo[1] if dm_lo else np.nan, "u2_max_dm_hi": dm_hi[1] if dm_hi else np.nan,
        "width_delta": 0.1, "had_frac": 0.7, "decay_samples": 50, "hit_estimator": "exact",
        "event_chunk": 1000, "seed_salt": np.nan,
    }


def test_member_is_the_rerun_shift_rebased_onto_the_band_central():
    band = pd.DataFrame([
        # legs named by width direction: at 2 GeV the "lo" leg sits above the "hi" leg
        _row("Ue", 1.0, 1e-8, 1e-4, dm_lo=(0.9e-8, 0.9e-4), dm_hi=(1.1e-8, 1.1e-4)),
        _row("Ue", 2.0, 2e-8, 2e-4, dm_lo=(2.2e-8, 2.2e-4), dm_hi=(1.8e-8, 1.8e-4)),
        _row("Ue", 3.0, 3e-8, 3e-4, dm_lo=(np.nan, np.nan), dm_hi=(2.7e-8, 2.7e-4)),
    ])
    control = pd.DataFrame([
        _row("Ue", 1.0, 1.001e-8, 1.001e-4),                     # same production, tiny code drift
        _row("Ue", 2.0, 2.5e-8, 2.5e-4, n_hits=101),             # different production
        _row("Ue", 3.0, 3e-8, 3e-4),
    ])
    exhad = pd.DataFrame([
        _row("Ue", 1.0, 1.001e-8 * 10**-0.3, 1.001e-4 * 10**0.2),
        _row("Ue", 2.0, 2.5e-8 * 10**0.05, 2.5e-4 * 10**-0.05),
        _row("Ue", 3.0, 3e-8 * 10**-0.5, 3e-4),
        _row("Ue", 3.7, 4e-8, 4e-4),                               # restored sensitivity
    ])
    out, topology = dmx.combine(band, control, exhad)
    r1, r2, r3 = (out[np.isclose(out.mass_GeV, m)].iloc[0] for m in (1.0, 2.0, 3.0))

    # the shift is control -> exHad, applied to the BAND central, not to the control
    assert r1["u2_min_exhad_shift_dex"] == pytest.approx(-0.3)
    assert r1["u2_min_exhad"] == pytest.approx(1e-8 * 10**-0.3)
    assert r1["u2_max_exhad"] == pytest.approx(1e-4 * 10**0.2)
    assert r1["u2_min_control_ratio"] == pytest.approx(1.001)
    assert bool(r1["same_production_as_band"])
    # the widened band is ordered and spans central, both legs and the member
    assert r1["u2_min_dm_lo"] == pytest.approx(1e-8 * 10**-0.3)
    assert r1["u2_min_dm_hi"] == pytest.approx(1.1e-8)
    assert r1["u2_max_dm_lo"] == pytest.approx(0.9e-4)
    assert r1["u2_max_dm_hi"] == pytest.approx(1e-4 * 10**0.2)
    # the published legs survive verbatim, even where they cross
    assert r2["u2_min_width_lo"] == pytest.approx(2.2e-8) and r2["u2_min_width_hi"] == pytest.approx(1.8e-8)
    assert r2["u2_min_exhad"] == pytest.approx(2e-8 * 10**0.05)          # 2.244e-8, above both legs
    assert r2["u2_min_dm_lo"] == pytest.approx(1.8e-8)
    assert r2["u2_min_dm_hi"] == pytest.approx(2e-8 * 10**0.05)
    assert r2["u2_max_dm_lo"] == pytest.approx(2e-4 * 10**-0.05) and r2["u2_max_dm_hi"] == pytest.approx(2.2e-4)
    assert not bool(r2["same_production_as_band"])
    # a leg without the boundary keeps the published convention: no ribbon
    assert np.isnan(r3["u2_min_dm_lo"]) and np.isnan(r3["u2_min_dm_hi"])
    assert r3["u2_min_exhad"] == pytest.approx(3e-8 * 10**-0.5)     # the member itself is still recorded
    assert not bool(r3["u2_min_exhad_member_missing"])
    assert topology["exhad_sensitive_where_central_is_not"] == [["Ue", 3.7, 4e-8, 4e-4]]
    assert topology["control_without_central_boundary"] == []
    assert list(out.columns[:len(band.columns)]) == list(band.columns)


def test_missing_member_is_flagged_not_invented():
    band = pd.DataFrame([_row("Umu", 3.62, 4e-7, 5e-7, dm_lo=(np.nan, np.nan), dm_hi=(1.8e-7, 1e-6))])
    control = pd.DataFrame([_row("Umu", 3.62, np.nan, np.nan, has_sensitivity=False)])
    exhad = pd.DataFrame([_row("Umu", 3.62, 2e-7, 7e-7)])
    out, topology = dmx.combine(band, control, exhad)
    row = out.iloc[0]
    assert bool(row["u2_min_exhad_member_missing"]) and np.isnan(row["u2_min_exhad"])
    assert np.isnan(row["u2_min_exhad_shift_dex"])
    assert topology["control_without_central_boundary"] == [["Umu", 3.62, "u2_min"], ["Umu", 3.62, "u2_max"]]
    assert row["n_hits_control"] == 100 and bool(row["same_production_as_band"])


def test_open_upper_edge_has_no_member():
    band = pd.DataFrame([_row("Ue", 0.5, 1e-7, np.nan, u2_max_open=True, dm_lo=(0.9e-7, np.nan), dm_hi=(1.1e-7, np.nan))])
    control = pd.DataFrame([_row("Ue", 0.5, 1e-7, np.nan, u2_max_open=True)])
    exhad = pd.DataFrame([_row("Ue", 0.5, 1.2e-7, np.nan, u2_max_open=True)])
    out, _ = dmx.combine(band, control, exhad)
    row = out.iloc[0]
    assert row["u2_min_exhad_shift_dex"] == pytest.approx(np.log10(1.2))
    assert np.isnan(row["u2_max_exhad"]) and np.isnan(row["u2_max_dm_lo"])
    assert not bool(row["u2_max_exhad_member_missing"])        # nothing to add on an open edge
    assert row["u2_min_dm_hi"] == pytest.approx(1.2e-7)
