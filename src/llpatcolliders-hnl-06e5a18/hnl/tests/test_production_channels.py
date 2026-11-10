"""Regression tests for the induced-tau pool and the b-baryon channel.

These exercise the production interfaces (REMAINING_WORK item 21): the
11-source induced-tau table with its polarization tags, the per-block
process_flavor path, and the Lambda_b driver through the dq2dm122 sampler.
The branching-ratio comparisons use a loose tolerance because HNLCalc's
integrator is Monte Carlo (nsample=500, ~5% error); they pin factor-2 and
|U|^2-vs-|U|^4 normalization bugs, not the rate model itself (item 13).
"""

import numpy as np

from production.constants import (
    M_TAU, M_LAMBDA_B, M_LAMBDA_C, FRAG_LAMBDA_B, LEPTON_MASSES,
)
from production.decay_engine import generate_baryon_csvs as gb
from production.decay_engine import generate_induced_tau as gi
from production.decay_engine.tau_decay import (
    TAU_MESON_2BODY_ASYMMETRY, compute_tau_production_br_components,
)
from production.hnlcalc import init_hnlcalc
from production.io import llp_csv_path


def test_tau_sources_inventory():
    """11 sources; leptonic 2-body ones carry the wrong-helicity asymmetry,
    semitauonic/baryonic ones are unpolarized; all kinematically open."""
    assert len(gi.TAU_SOURCES) == 11
    kinds = [s[4] for s in gi.TAU_SOURCES]
    assert kinds.count("2body") == 4
    assert kinds.count("baryon") == 1
    for parent, m_parent, dau, mdau, kind, quark, frag, br, asym in gi.TAU_SOURCES:
        assert 0.0 < br < 0.2, f"{parent}: BR={br}"
        assert 0.0 < frag <= 1.0, f"{parent}: frag={frag}"
        if kind == "2body":
            assert m_parent > M_TAU, f"{parent}: tau nu closed"
            assert asym == TAU_MESON_2BODY_ASYMMETRY
        else:
            assert m_parent > mdau + M_TAU, f"{parent} -> {dau} tau nu closed"
            assert asym == 0.0


def test_induced_tau_process_flavor(tmp_path, monkeypatch):
    """Synthetic tau pool through process_flavor: row count preserved across
    polarization blocks, weight = tau_w * BR(tau -> N X), sampled HNLs on
    shell, zero-byte sentinel above m_tau."""
    monkeypatch.setattr(gi, "OUTPUT_BASE", tmp_path)
    rng = np.random.default_rng(7)
    n = 400
    p = rng.uniform(5.0, 50.0, n)
    cos_t = rng.uniform(0.9, 1.0, n)
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    sin_t = np.sqrt(1.0 - cos_t**2)
    tau_4v = np.column_stack([
        np.sqrt(p**2 + M_TAU**2),
        p * sin_t * np.cos(phi), p * sin_t * np.sin(phi), p * cos_t,
    ])
    tau_w = np.full(n, 2.5e3)
    tau_asym = np.concatenate([np.full(n // 2, -1.0), np.zeros(n - n // 2)])
    m_N = 0.5

    gi.process_flavor("Umu", tau_4v, tau_w, tau_asym, [m_N, 2.0], rng)

    assert llp_csv_path("Umu", "induced_tau", 2.0, base=tmp_path).stat().st_size == 0

    data = np.loadtxt(llp_csv_path("Umu", "induced_tau", m_N, base=tmp_path),
                      delimiter=",")
    assert data.shape == (n, 5)

    # One common weight tau_w * br_total across both polarization blocks.
    weights = np.unique(data[:, 0])
    assert len(weights) == 1 and weights[0] > 0
    br_implied = weights[0] / tau_w[0]
    assert 0.0 < br_implied < 1.0
    _, _, br_ref = compute_tau_production_br_components(init_hnlcalc("Umu"), m_N)
    assert np.isclose(br_implied, br_ref, rtol=0.3)

    # Sampled HNLs are on shell at m_N.
    m2 = data[:, 1]**2 - data[:, 2]**2 - data[:, 3]**2 - data[:, 4]**2
    assert np.allclose(m2, m_N**2, atol=1e-3)


def test_baryon_channel_through_production_interface(tmp_path, monkeypatch):
    """Lambda_b -> Lambda_c l N via decay_3body_weighted_dq2dm122: on-shell
    HNLs, weight normalization 2*sigma_b*f_Lb*BR, sentinel above threshold."""
    monkeypatch.setattr(gb, "OUTPUT_BASE", tmp_path)
    rng = np.random.default_rng(11)
    n = 300
    pool, sigma_b = gb.build_lambda_b_pool(n, rng)
    m_N = 1.0

    gb.process_flavor("Ue", pool, sigma_b, [m_N, 4.0], rng)

    # 4.0 GeV > M_LAMBDA_B - M_LAMBDA_C - m_e: kinematically closed sentinel.
    assert llp_csv_path("Ue", "Bbaryon", 4.0, base=tmp_path).stat().st_size == 0

    data = np.loadtxt(llp_csv_path("Ue", "Bbaryon", m_N, base=tmp_path),
                      delimiter=",")
    assert data.shape == (n, 5)
    assert (data[:, 0] > 0).all()

    m2 = data[:, 1]**2 - data[:, 2]**2 - data[:, 3]**2 - data[:, 4]**2
    assert np.allclose(m2, m_N**2, atol=1e-3)

    hnl = init_hnlcalc("Ue")
    dbr = hnl.get_3body_dbr_baryon("5122", "4122", "11")
    br_ref = hnl.integrate_3body_br(
        dbr, m_N, M_LAMBDA_B, M_LAMBDA_C, LEPTON_MASSES["Ue"],
        coupling=1.0, nsample=500, integration="dq2dm122",
    )
    assert np.isclose(data[:, 0].sum(), 2.0 * sigma_b * FRAG_LAMBDA_B * br_ref,
                      rtol=0.3)
