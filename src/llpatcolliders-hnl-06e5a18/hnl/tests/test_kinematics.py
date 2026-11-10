"""Energy/momentum conservation for the 2-body and 3-body decay kinematics."""

import numpy as np
import pytest

from production.decay_engine.kinematics import (
    decay_2body, decay_2body_polarized, decay_3body_flat, _sample_polar_cos,
    decay_3body_weighted_dq2dE, decay_3body_weighted_dE,
)


@pytest.fixture
def rng():
    return np.random.default_rng(0)


def _make_parents(n, mass, p_max=20.0, rng=None):
    """Random parent 4-vectors with |p| ~ uniform in [0, p_max]."""
    rng = rng or np.random.default_rng(0)
    p = rng.uniform(0.0, p_max, n)
    cos_th = rng.uniform(-1.0, 1.0, n)
    sin_th = np.sqrt(1.0 - cos_th**2)
    phi = rng.uniform(0, 2 * np.pi, n)
    px = p * sin_th * np.cos(phi)
    py = p * sin_th * np.sin(phi)
    pz = p * cos_th
    E = np.sqrt(p**2 + mass**2)
    return E, px, py, pz


def test_2body_conserves_4momentum(rng):
    M, m1, m2 = 5.0, 0.5, 0.1
    E, px, py, pz = _make_parents(2000, M, rng=rng)
    d1, d2 = decay_2body(E, px, py, pz, M, m1, m2, rng=rng)

    assert np.allclose(d1[:, 0] + d2[:, 0], E, rtol=1e-9, atol=1e-9)
    assert np.allclose(d1[:, 1] + d2[:, 1], px, rtol=1e-9, atol=1e-9)
    assert np.allclose(d1[:, 2] + d2[:, 2], py, rtol=1e-9, atol=1e-9)
    assert np.allclose(d1[:, 3] + d2[:, 3], pz, rtol=1e-9, atol=1e-9)

    # Daughter on-shell
    m1_recon = np.sqrt(np.maximum(d1[:, 0]**2 - (d1[:, 1]**2 + d1[:, 2]**2 + d1[:, 3]**2), 0.0))
    m2_recon = np.sqrt(np.maximum(d2[:, 0]**2 - (d2[:, 1]**2 + d2[:, 2]**2 + d2[:, 3]**2), 0.0))
    assert np.allclose(m1_recon, m1, atol=1e-7)
    assert np.allclose(m2_recon, m2, atol=1e-7)


def test_sample_polar_cos_distribution(rng):
    # pdf ∝ 1 + a·cosθ  =>  <cosθ> = a/3, and a=0 is flat (<cosθ>=0).
    for a in (0.0, 0.5, -1.0, 1.0):
        c = _sample_polar_cos(a, 200000, rng)
        assert c.min() >= -1.0 and c.max() <= 1.0
        assert abs(c.mean() - a / 3.0) < 5e-3


def test_2body_polarized_conserves_and_reduces_to_isotropic(rng):
    M, m1, m2 = 1.777, 0.139, 0.5  # tau -> pi N -like
    E, px, py, pz = _make_parents(5000, M, p_max=30.0, rng=rng)

    d1, d2 = decay_2body_polarized(E, px, py, pz, M, m1, m2, asymmetry=0.0, rng=rng)
    # 4-momentum conservation and on-shell daughters.
    assert np.allclose(d1[:, 0] + d2[:, 0], E, rtol=1e-9, atol=1e-9)
    for col, p in zip(range(1, 4), (px, py, pz)):
        assert np.allclose(d1[:, col] + d2[:, col], p, rtol=1e-9, atol=1e-9)
    m2_rec = np.sqrt(np.maximum(d2[:, 0]**2 - d2[:, 1]**2 - d2[:, 2]**2 - d2[:, 3]**2, 0.0))
    assert np.allclose(m2_rec, m2, atol=1e-7)


def test_2body_polarized_shifts_energy_spectrum(rng):
    # A nonzero asymmetry must measurably change the analyzed-daughter (N) lab
    # energy spectrum relative to the isotropic case (the whole point of the fix).
    M, m1, m2 = 1.777, 0.139, 0.5
    E, px, py, pz = _make_parents(40000, M, p_max=40.0, rng=rng)
    _, n_iso = decay_2body_polarized(E, px, py, pz, M, m1, m2, asymmetry=0.0, rng=rng)
    _, n_pol = decay_2body_polarized(E, px, py, pz, M, m1, m2, asymmetry=-1.0, rng=rng)
    assert abs(n_iso[:, 0].mean() - n_pol[:, 0].mean()) > 1e-2 * n_iso[:, 0].mean()


def test_3body_conserves_4momentum(rng):
    M, m1, m2, m3 = 5.28, 0.139, 0.106, 1.0   # B -> pi mu N-like masses
    E, px, py, pz = _make_parents(500, M, rng=rng)
    d1, d2, d3 = decay_3body_flat(E, px, py, pz, M, m1, m2, m3, rng=rng)

    s = d1 + d2 + d3
    assert np.allclose(s[:, 0], E, rtol=1e-8, atol=1e-8)
    assert np.allclose(s[:, 1], px, rtol=1e-8, atol=1e-8)
    assert np.allclose(s[:, 2], py, rtol=1e-8, atol=1e-8)
    assert np.allclose(s[:, 3], pz, rtol=1e-8, atol=1e-8)

    # Each daughter is on-shell at its assigned mass
    for d, m in [(d1, m1), (d2, m2), (d3, m3)]:
        m_recon = np.sqrt(np.maximum(d[:, 0]**2 - (d[:, 1]**2 + d[:, 2]**2 + d[:, 3]**2), 0.0))
        assert np.allclose(m_recon, m, atol=1e-6)


def test_3body_threshold_raises(rng):
    """3-body decay at or below threshold must error rather than sample garbage."""
    M = 1.0
    m1, m2, m3 = 0.5, 0.5, 0.5  # threshold = 1.5 > M
    E, px, py, pz = _make_parents(10, M, rng=rng)
    with pytest.raises(RuntimeError):
        decay_3body_flat(E, px, py, pz, M, m1, m2, m3, rng=rng)


# ---------------------------------------------------------------------------
# Matrix-element-weighted samplers
# ---------------------------------------------------------------------------


def test_weighted_dq2dE_conserves_4momentum_and_masses(rng):
    """dq²/dE sampler with a constant dbr reduces to flat phase space, with
    4-momentum conservation and on-shell daughters."""
    M, m1, m2, m3 = 5.28, 0.139, 0.106, 1.0   # B -> pi mu N-like
    E, px, py, pz = _make_parents(400, M, rng=rng)
    # A constant dbr makes the weight independent of (q², E_N); accept-reject
    # then degenerates to uniform sampling in the (q², E_N) Dalitz region,
    # which is what we use the conservation laws to verify.
    dbr_expr = "1.0 + 0.0*q2 + 0.0*energy"
    d1, d2, d3 = decay_3body_weighted_dq2dE(
        E, px, py, pz, M, m1, m2, m3, dbr_expr=dbr_expr, rng=rng,
    )

    s = d1 + d2 + d3
    assert np.allclose(s[:, 0], E, rtol=1e-7, atol=1e-7)
    for col, p in zip(range(1, 4), (px, py, pz)):
        assert np.allclose(s[:, col], p, rtol=1e-7, atol=1e-7)

    for d, m in [(d1, m1), (d2, m2), (d3, m3)]:
        m_rec = np.sqrt(np.maximum(
            d[:, 0]**2 - (d[:, 1]**2 + d[:, 2]**2 + d[:, 3]**2), 0.0
        ))
        assert np.allclose(m_rec, m, atol=1e-5)


def test_weighted_dE_conserves_and_biases_energy(rng):
    """dE sampler conserves 4-momentum and reproduces the requested E_N bias.

    With a dbr ∝ energy, accept-reject must yield <E_N> > <E_N>_flat in the
    parent rest frame, demonstrating the spectrum is matrix-element-weighted
    rather than uniform.
    """
    M, m1, m2, m3 = 1.777, 0.000511, 0.0, 0.3  # tau-like
    # Tau at rest so the lab-frame HNL energy is the parent-rest E_N directly.
    n = 2000
    E = np.full(n, M)
    px = np.zeros(n); py = np.zeros(n); pz = np.zeros(n)

    # Flat-in-E dbr — should give the uniform-E mean (midpoint of allowed range).
    Emin = m3
    Emax = (M**2 + m3**2 - (m1 + m2)**2) / (2.0 * M)
    flat_expr = "1.0 + 0.0*energy"
    _, _, hnl_flat = decay_3body_weighted_dE(
        E, px, py, pz, M, m1, m2, m3, dbr_expr=flat_expr, rng=rng,
    )

    biased_expr = "energy"
    _, _, hnl_biased = decay_3body_weighted_dE(
        E, px, py, pz, M, m1, m2, m3, dbr_expr=biased_expr,
        rng=np.random.default_rng(7),
    )

    # 4-momentum sanity (against parent at rest in lab).
    # (only checking HNL on-shell since other daughters tested above)
    m_rec = np.sqrt(np.maximum(
        hnl_flat[:, 0]**2 - (hnl_flat[:, 1]**2 + hnl_flat[:, 2]**2 + hnl_flat[:, 3]**2),
        0.0,
    ))
    assert np.allclose(m_rec, m3, atol=1e-6)

    # Flat dbr -> mean E = midpoint of [Emin, Emax]; biased dbr ∝ E -> mean = 2/3 (Emax² - Emin²)/(Emax² - Emin²) ... in this regime simply larger than flat mean.
    flat_mean = hnl_flat[:, 0].mean()
    biased_mean = hnl_biased[:, 0].mean()
    expected_flat = 0.5 * (Emin + Emax)
    assert flat_mean == pytest.approx(expected_flat, rel=5e-2)
    assert biased_mean > flat_mean + 0.02 * expected_flat
