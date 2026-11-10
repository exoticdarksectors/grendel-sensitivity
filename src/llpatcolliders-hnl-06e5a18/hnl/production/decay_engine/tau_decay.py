"""Shared tau -> N + X branching fractions and decay sampler."""

import numpy as np

from production.constants import M_TAU
from production.decay_engine.kinematics import (
    decay_2body_polarized, decay_3body_weighted_dE,
)
from production.hnlcalc import init_hnlcalc


# Longitudinal asymmetry of the N direction relative to the tau lab momentum
# in the 2-body tau -> meson N decay (pdf ~ 1 + asym*cos(theta_N)), in the
# fixed unit-analyzing-power approximation. The sign is set by the tau
# helicity at production and is the same for both tau charges within a source:
#  - W -> tau nu yields natural-helicity taus (h = -1 for tau-, +1 for tau+);
#    the N carries the spin projection and is emitted forward -> asym = +1.
#  - P+ -> tau+ nu leptonic decays (Ds/D/B/Bc) yield helicity-suppressed
#    "wrong-helicity" taus (h = -1 for tau+, +1 for tau-); the N is emitted
#    backward -> asym = -1.
TAU_W_2BODY_ASYMMETRY = +1.0
TAU_MESON_2BODY_ASYMMETRY = -1.0


TAU_2BODY_MESON_PDGS = [211, 321, 213, 323]   # pi+, K+, rho+, K*+ (HNLCalc convention)


def _eval_tau_2body_br(hnl, meson_pdg, m_N):
    br_expr = hnl.get_2body_br_tau(15, meson_pdg)
    br_val = eval(br_expr, {"np": np, "__builtins__": {}},
                  {"mass": m_N, "coupling": 1.0})
    if np.isnan(br_val) or br_val < 0:
        return 0.0
    return float(br_val)


def _eval_tau_3body_br(hnl, lep_pid, nu_pid, m_N):
    dbr = hnl.get_3body_dbr_tau(15, -lep_pid, nu_pid)
    m_lep = hnl.masses(lep_pid)
    br_val = hnl.integrate_3body_br(
        dbr, m_N, M_TAU, m_lep, 0.0,
        coupling=1.0, nsample=500, integration="dE",
    )
    if br_val is None or np.isnan(br_val) or br_val < 0:
        return 0.0, None
    return float(br_val), dbr


def compute_tau_production_br_components(hnl, m_N):
    br_2body_channels = []
    for meson_pdg in TAU_2BODY_MESON_PDGS:
        if m_N >= M_TAU - hnl.masses(meson_pdg):
            continue
        br = _eval_tau_2body_br(hnl, meson_pdg, m_N)
        if br > 0:
            br_2body_channels.append((hnl.masses(meson_pdg), br))

    br_3body_channels = []
    for lep_pid in [11, 13]:
        m_lep = hnl.masses(lep_pid)
        if m_N >= M_TAU - m_lep:
            continue
        br_nt, dbr_nt = _eval_tau_3body_br(hnl, lep_pid, 16, m_N)
        if br_nt > 0 and dbr_nt is not None:
            br_3body_channels.append((m_lep, dbr_nt, br_nt))
        br_nl, dbr_nl = _eval_tau_3body_br(hnl, lep_pid, lep_pid + 1, m_N)
        if br_nl > 0 and dbr_nl is not None:
            br_3body_channels.append((m_lep, dbr_nl, br_nl))

    br_total = (
        sum(br for _, br in br_2body_channels)
        + sum(ch[2] for ch in br_3body_channels)
    )
    return br_2body_channels, br_3body_channels, br_total


def sample_hnl_from_tau(tau_E, tau_px, tau_py, tau_pz, m_N,
                        br_2body_channels, br_3body_channels, br_total, rng,
                        asymmetry):
    n_events = len(tau_E)
    hnl_4v = np.empty((n_events, 4))

    br_2body_total = sum(br for _, br in br_2body_channels)
    p_2body = br_2body_total / br_total if br_total > 0 else 0.0
    use_2body = rng.random(n_events) < p_2body
    use_3body = ~use_2body

    if use_2body.any() and br_2body_channels:
        br_arr = np.array([br for _, br in br_2body_channels], dtype=float)
        prob = br_arr / br_arr.sum()
        ch_idx = rng.choice(len(br_2body_channels), size=use_2body.sum(), p=prob)
        evt_idx = np.where(use_2body)[0]
        for k in np.unique(ch_idx):
            m_meson, _ = br_2body_channels[int(k)]
            sel = evt_idx[ch_idx == k]
            _, hnl_2b = decay_2body_polarized(
                tau_E[sel], tau_px[sel], tau_py[sel], tau_pz[sel],
                M_TAU, m_meson, m_N, asymmetry=asymmetry, rng=rng,
            )
            hnl_4v[sel] = hnl_2b

    if use_3body.any() and br_3body_channels:
        br_arr = np.array([ch[2] for ch in br_3body_channels], dtype=float)
        prob = br_arr / br_arr.sum()
        ch_idx = rng.choice(len(br_3body_channels), size=use_3body.sum(), p=prob)
        evt_idx = np.where(use_3body)[0]
        for k in np.unique(ch_idx):
            m_lep, dbr_expr, _ = br_3body_channels[int(k)]
            sel = evt_idx[ch_idx == k]
            _, _, hnl_3b = decay_3body_weighted_dE(
                tau_E[sel], tau_px[sel], tau_py[sel], tau_pz[sel],
                M_TAU, m_lep, 0.0, m_N,
                dbr_expr=dbr_expr, coupling=1.0, rng=rng,
            )
            hnl_4v[sel] = hnl_3b

    return hnl_4v, use_2body, use_3body
