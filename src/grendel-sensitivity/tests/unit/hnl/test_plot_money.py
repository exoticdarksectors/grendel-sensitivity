"""Regression tests for topology-safe uncertainty-ribbon interpolation."""

import numpy as np

from grendel.models.hnl.plot_money import _dex_densify, _metadata, _provenance


def test_dex_densify_does_not_extrapolate_past_finite_support():
    bm = np.array([1.0, 2.0, 3.0])
    central = np.ones(3)
    lower = np.array([0.8, 0.8, np.nan])
    upper = np.array([1.2, 1.2, np.nan])
    cm = np.array([1.0, 1.5, 2.0, 2.5, 3.0])

    _, lo, hi = _dex_densify(
        bm, central, lower, upper, cm, np.ones(len(cm)))

    assert np.all(np.isfinite(lo[:3]))
    assert np.all(np.isfinite(hi[:3]))
    assert np.all(np.isnan(lo[3:]))
    assert np.all(np.isnan(hi[3:]))
    assert np.all(lo[:3] < 1.0)
    assert np.all(hi[:3] > 1.0)


def test_dex_densify_does_not_bridge_invalid_topology_anchor():
    bm = np.arange(1.0, 6.0)
    central = np.ones(5)
    lower = np.array([0.8, 0.8, np.nan, 0.9, 0.9])
    upper = np.array([1.2, 1.2, np.nan, 1.1, 1.1])
    cm = np.arange(1.0, 5.1, 0.5)

    _, lo, hi = _dex_densify(
        bm, central, lower, upper, cm, np.ones(len(cm)))

    gap = (cm > 2.0) & (cm < 4.0)
    assert np.all(np.isnan(lo[gap]))
    assert np.all(np.isnan(hi[gap]))
    assert np.all(np.isfinite(lo[~gap]))
    assert np.all(np.isfinite(hi[~gap]))


def test_metadata_and_provenance_record_kaon_transport_band():
    metadata = _metadata("test", 3000.0, 100, have_fonll=True, have_kaon=True)

    assert any("charged-kaon transport" in item
               for item in metadata["scope"]["in_band"])
    assert "kaon_lower_edge" in metadata["band_sources"]
    assert not any("kaon flux/transport" in item
                   for item in metadata["limitations_not_banded"])

    provenance = _provenance(have_fonll=True, have_kaon=True)
    assert "charged-kaon transport" in provenance
    assert "purple" in provenance
