"""Unit tests for the FairShip->GRENDEL signal acceptance core (Stage 2)."""
import sys
from pathlib import Path

import numpy as np
import pytest

from grendel.reco import acceptance as dra
from grendel.geometry.grendel_geometry import mesh_fiducial  # noqa: E402


def test_boost_rest_to_lab_at_rest_is_identity():
    px = np.array([0.3, -0.3]); py = np.array([0.1, -0.1])
    pz = np.array([0.2, -0.2]); E = np.sqrt(px**2 + py**2 + pz**2 + 0.01)
    lx, ly, lz, lE = dra.boost_rest_to_lab([1.0, 0, 0, 0], px, py, pz, E)
    assert np.allclose(lx, px) and np.allclose(ly, py)
    assert np.allclose(lz, pz) and np.allclose(lE, E)


def test_boost_rest_to_lab_increases_energy_and_conserves_momentum():
    # consistent 2-body decay of a mass-1 parent: back-to-back, each E = M/2,
    # so the rest-frame daughter energies sum to the parent mass.
    px = np.array([0.5, -0.5]); py = np.array([0.0, 0.0]); pz = np.array([0.0, 0.0])
    E = np.array([0.5, 0.5])
    p4 = np.array([10.0, 0.0, 0.0, np.sqrt(100 - 1.0)])  # parent (m=1) moving in +z
    lx, ly, lz, lE = dra.boost_rest_to_lab(p4, px, py, pz, E)
    assert lE.sum() == pytest.approx(p4[0], rel=1e-9)          # energy conserved
    assert lz.sum() == pytest.approx(p4[3], rel=1e-9)          # pz conserved
    assert lE.sum() > E.sum()                                  # boosted up


def test_best_two_directions_picks_leading_charged():
    # 3 charged + 1 neutral; expect the two highest-|p| charged, unit dirs
    px = np.array([5.0, 0.0, 3.0, 0.0]); py = np.array([0.0, 4.0, 0.0, 9.0])
    pz = np.array([0.0, 0.0, 0.0, 0.0])
    E = np.sqrt(px**2 + py**2 + pz**2)      # massless -> beta = 1
    charge = np.array([1.0, -1.0, 1.0, 0.0]); stable = np.array([1, 1, 1, 1])
    out = dra.best_two_directions(px, py, pz, E, charge, stable)
    assert out is not None
    d1, d2, p_soft, b1, b2 = out
    assert np.allclose(np.linalg.norm(d1), 1.0)
    assert np.allclose(d1, [1, 0, 0])          # |p|=5 leading
    assert np.allclose(d2, [0, 1, 0])          # |p|=4 second (neutral 9 excluded)
    assert p_soft == pytest.approx(4.0)
    assert b1 == pytest.approx(1.0) and b2 == pytest.approx(1.0)


def test_best_two_directions_returns_true_beta_for_massive_tracks():
    # two charged pions at |p| = 0.2 GeV: beta = p/E = 0.2/sqrt(0.2^2+m_pi^2)
    m_pi = 0.13957
    px = np.array([0.2, -0.2]); py = np.zeros(2); pz = np.zeros(2)
    E = np.sqrt(px**2 + m_pi**2)
    charge = np.array([1.0, -1.0]); stable = np.array([1, 1])
    out = dra.best_two_directions(px, py, pz, E, charge, stable)
    assert out is not None
    _, _, _, b1, b2 = out
    beta_true = 0.2 / np.sqrt(0.2**2 + m_pi**2)   # ~0.82
    assert b1 == pytest.approx(beta_true, rel=1e-9)
    assert b2 == pytest.approx(beta_true, rel=1e-9)


def test_best_two_directions_needs_two_above_cut():
    px = np.array([5.0, 0.1]); py = np.zeros(2); pz = np.zeros(2)
    E = np.sqrt(px**2 + py**2 + pz**2)
    charge = np.array([1.0, -1.0]); stable = np.array([1, 1])
    assert dra.best_two_directions(px, py, pz, E, charge, stable) is None  # 2nd below P_CUT


def _hitting_vertex(rng):
    origin = np.zeros(3)
    for _ in range(5000):
        d = rng.normal(size=3); d[1] = abs(d[1]) + 0.5; d /= np.linalg.norm(d)
        loc, _, _ = mesh_fiducial.ray.intersects_location([origin], [d])
        if len(loc) >= 2:
            ds = np.sort(np.linalg.norm(loc - origin, axis=1))
            return origin + 0.5 * (ds[0] + ds[1]) * d, d
    pytest.skip("no hitting direction found")


def test_reconstruct_and_select_smoke():
    rng = np.random.default_rng(1)
    vtx, d = _hitting_vertex(rng)
    # two forward daughters with a small opening about the HNL direction
    e1 = np.cross(d, [0, 0, 1.0]); e1 /= np.linalg.norm(e1)
    dir1 = d + 0.05 * e1; dir1 /= np.linalg.norm(dir1)
    dir2 = d - 0.05 * e1; dir2 /= np.linalg.norm(dir2)
    N = 16
    mc = dra.reconstruct_decays(
        np.tile(vtx, (N, 1)), np.tile(dir1, (N, 1)), np.tile(dir2, (N, 1)),
        np.full(N, 5.0), dra.HIT_RESOLUTION, dra.SIGMA_T_DEFAULT, rng)
    for k in ('sep', 'sep_outer', 'open_angle', 'dca', 'collin', 'pointing',
              'vtx_in', 'on_tracker', 'timing_chi2', 'p_soft'):
        assert k in mc and len(mc[k]) == N
    sel = dra.selection_mask(mc)
    assert sel.dtype == bool and len(sel) == N


def test_timing_chi2_penalizes_slow_daughters():
    # Same geometry, same smear noise (fresh identically-seeded rng per call):
    # daughters at beta = 0.5 arrive late relative to the c-hypothesis, so the
    # timing chi2 must grow vs the beta = 1 baseline. Guards the fix that
    # threads real daughter betas into timing_chi2_4hit.
    rng0 = np.random.default_rng(3)
    vtx, d = _hitting_vertex(rng0)
    e1 = np.cross(d, [0, 0, 1.0]); e1 /= np.linalg.norm(e1)
    dir1 = d + 0.05 * e1; dir1 /= np.linalg.norm(dir1)
    dir2 = d - 0.05 * e1; dir2 /= np.linalg.norm(dir2)
    N = 32
    args = (np.tile(vtx, (N, 1)), np.tile(dir1, (N, 1)), np.tile(dir2, (N, 1)),
            np.full(N, 5.0), dra.HIT_RESOLUTION, dra.SIGMA_T_DEFAULT)
    mc_fast = dra.reconstruct_decays(*args, np.random.default_rng(7))
    mc_slow = dra.reconstruct_decays(*args, np.random.default_rng(7),
                                     beta1=np.full(N, 0.5),
                                     beta2=np.full(N, 0.5))
    ok = np.isfinite(mc_fast['timing_chi2']) & np.isfinite(mc_slow['timing_chi2'])
    if not ok.any():
        pytest.skip("no reconstructable decays for this geometry draw")
    assert (np.median(mc_slow['timing_chi2'][ok])
            > np.median(mc_fast['timing_chi2'][ok]))


def _synthetic_templates():
    """Two decay templates, each e+ e- back-to-back (2 charged stable tracks)."""
    # rest-frame: each daughter |p| = 0.5 GeV (m_N ~ 1 GeV), opposite directions
    pdg = np.array([11, -11, 11, -11], dtype=np.int32)
    px = np.array([0.5, -0.5, 0.0, 0.0]); py = np.array([0.0, 0.0, 0.5, -0.5])
    pz = np.zeros(4); energy = np.full(4, 0.5); mass = np.zeros(4)
    charge = np.array([-1.0, 1.0, -1.0, 1.0]); stable = np.ones(4, bool)
    return dict(daughter_counts=np.array([2, 2], np.int32),
                pdg=pdg, px=px, py=py, pz=pz, energy=energy, mass=mass,
                charge=charge, stable=stable)


def test_scan_u2_has_interior_lifetime_peak():
    rng = np.random.default_rng(2)
    origin = np.zeros(3)
    p4s, dirs, ent, exi = [], [], [], []
    m = 1.0
    while len(p4s) < 8:
        d = rng.normal(size=3); d[1] = abs(d[1]) + 0.4; d /= np.linalg.norm(d)
        loc, _, _ = mesh_fiducial.ray.intersects_location([origin], [d])
        if len(loc) < 2:
            continue
        ds = np.sort(np.linalg.norm(loc - origin, axis=1))
        P = 30.0; E = np.sqrt(P * P + m * m)
        p4s.append([E, *(P * d)]); dirs.append(d); ent.append(ds[0]); exi.append(ds[1])
    p4s = np.array(p4s); dirs = np.array(dirs); ent = np.array(ent); exi = np.array(exi)
    bg = np.array([np.linalg.norm(p[1:]) / m for p in p4s])
    T = _synthetic_templates()
    d, passed, _tmpl_idx = dra.build_event_mc(p4s, dirs, ent, exi, T, n_samples=40, rng=rng)
    u2 = np.logspace(-9, -1, 30)
    u2, N = dra.scan_u2(d, passed, exi - ent, np.ones(len(p4s)), bg,
                        ctau_u2_1=1e-3, L_int_pb=3e6, u2_grid=u2)
    ipk = int(np.argmax(N))
    assert 0 < ipk < len(u2) - 1          # peak is interior (lifetime frontier)
    assert N[0] < N[ipk] and N[-1] < N[ipk]
    assert N[ipk] > 0


def test_signal_diagnostics_separate_sample_and_event_ess():
    d = np.array([[1.0, 1.5], [1.0, 1.5]])
    passed = np.array([[True, True], [True, True]])
    diagnostics = dra.signal_contribution_diagnostics(
        d,
        passed,
        path_len=np.array([1.0, 1.0]),
        weight=np.array([100.0, 1.0]),
        beta_gamma=np.array([1.0, 1.0]),
        ctau_u2_1=1.0,
        u2=1.0,
    )

    assert diagnostics["sample_ess"] > diagnostics["event_ess"]
    assert diagnostics["event_ess"] < 1.1
    assert diagnostics["max_event_fraction"] > 0.99
    assert diagnostics["nonzero_samples"] == 4
    assert diagnostics["nonzero_events"] == 2
