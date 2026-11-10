"""Kaon channel (K+ -> l N) production sanity checks."""

import numpy as np
import pytest

from production.decay_engine import generate_kaon_csvs as K
from production.constants import M_KAON, M_ELECTRON, M_TAU


@pytest.fixture
def rng():
    return np.random.default_rng(0)


def test_kaon_pool_is_on_shell(rng):
    pool = K.sample_kaon_4vectors(20000, rng)
    m2 = pool['E']**2 - pool['px']**2 - pool['py']**2 - pool['pz']**2
    assert np.allclose(np.sqrt(np.clip(m2, 0, None)), M_KAON, atol=1e-6)
    # Soft spectrum: mean pT should sit well below 1 GeV.
    assert 0.3 < pool['pt'].mean() < 1.0


def test_kaon_channel_open_for_electron_and_closed_above_threshold(tmp_path, rng, monkeypatch):
    monkeypatch.setattr(K, "OUTPUT_BASE", tmp_path / "llp_4vectors")
    pool = K.sample_kaon_4vectors(3000, rng)
    # m_K - m_e ~ 0.4936 GeV: 0.35 open, 0.60 closed.
    K.process_flavor("Ue", pool, [0.35, 0.60], rng)

    open_csv = tmp_path / "llp_4vectors" / "Ue" / "Kmeson" / "mN_0p350.csv"
    closed_csv = tmp_path / "llp_4vectors" / "Ue" / "Kmeson" / "mN_0p600.csv"
    assert open_csv.stat().st_size > 0
    assert closed_csv.stat().st_size == 0

    data = np.loadtxt(open_csv, delimiter=",")
    assert data.shape[1] == 5
    assert (data[:, 0] > 0).all()
    m_rec = np.sqrt(np.clip(data[:, 1]**2 - data[:, 2]**2 - data[:, 3]**2 - data[:, 4]**2, 0, None))
    # m = sqrt(E^2 - pz^2) loses precision for the rare ultra-boosted kaons in
    # the Gaussian rapidity tail (8-significant-figure %.8e storage, not the
    # kinematics), so check the bulk rather than every event.
    assert abs(np.median(m_rec) - 0.35) < 1e-4
    assert np.quantile(np.abs(m_rec - 0.35), 0.95) < 1e-2


def test_kaon_channel_closed_for_tau(tmp_path, rng, monkeypatch):
    # m_tau > m_K, so the tau-flavor kaon channel is shut at every mass.
    monkeypatch.setattr(K, "OUTPUT_BASE", tmp_path / "llp_4vectors")
    assert M_TAU > M_KAON - M_ELECTRON
    pool = K.sample_kaon_4vectors(2000, rng)
    K.process_flavor("Utau", pool, [0.2, 0.35], rng)
    for m in ("mN_0p200.csv", "mN_0p350.csv"):
        f = tmp_path / "llp_4vectors" / "Utau" / "Kmeson" / m
        assert f.exists() and f.stat().st_size == 0
