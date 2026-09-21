"""Regression tests for exact-hit, resumable HNL analysis campaigns."""

import numpy as np

from grendel.models.hnl.spec import seed_for as _seed_for, select_hit_sample as _select_hit_sample
from grendel.reco.acceptance import scan_u2
from grendel.band.hnl.campaign import run_tag as band_run_tag


def test_point_seed_is_stable_and_salt_creates_independent_stream():
    assert _seed_for("Ue", "1p000") == _seed_for("Ue", "1p000")
    assert _seed_for("Ue", "1p000") != _seed_for(
        "Ue", "1p000", "control-1")
    assert _seed_for("Ue", "1p000", "control-1") == _seed_for(
        "Ue", "1p000", "control-1")


def test_exact_hit_selection_retains_indices_and_weights():
    idx = np.array([1, 4, 8, 12])
    weights = np.array([0.1, 0.2, 0.3, 0.4])
    got_idx, got_weights, estimator = _select_hit_sample(
        idx, weights, None, np.random.default_rng(3))

    assert estimator == "exact"
    assert np.array_equal(got_idx, idx)
    assert np.array_equal(got_weights, weights)


def test_scan_u2_is_additive_over_event_chunks():
    rng = np.random.default_rng(4)
    n_events, n_samples = 9, 7
    d = rng.uniform(20.0, 30.0, size=(n_events, n_samples))
    passed = rng.random((n_events, n_samples)) > 0.35
    path_len = rng.uniform(3.0, 8.0, size=n_events)
    weight = rng.uniform(0.1, 2.0, size=n_events)
    beta_gamma = rng.uniform(2.0, 50.0, size=n_events)
    u2 = np.logspace(-9, -3, 13)
    kwargs = dict(ctau_u2_1=2.5e-4, L_int_pb=3e6, u2_grid=u2)

    _, full = scan_u2(
        d, passed, path_len, weight, beta_gamma, **kwargs)
    split = np.zeros_like(full)
    for start, stop in ((0, 4), (4, 7), (7, 9)):
        _, part = scan_u2(
            d[start:stop], passed[start:stop], path_len[start:stop],
            weight[start:stop], beta_gamma[start:stop], **kwargs)
        split += part

    assert np.allclose(split, full, rtol=2e-15, atol=0.0)


def test_band_run_tag_preserves_legacy_and_hashes_campaign_config():
    variation = {
        "name": "pdf_0017",
        "bottom": {"sha256": "a"},
        "charm": {"sha256": "b"},
    }
    assert band_run_tag(variation) == "band_pdf_0017_fb8e20fc"
    assert band_run_tag(variation, {"x": 1}) == (
        "band_pdf_0017_fb8e20fc_cfg5041bf1f")
    assert band_run_tag(variation, {"x": 1}) != band_run_tag(
        variation, {"x": 2})
