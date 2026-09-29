"""Scalar-portal model layer (PBC benchmark BC4)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from particle import Particle

G_F = 1.1663787e-5
HBAR = 6.582119569e-25
C_LIGHT = 299792458.0
HBAR_C = HBAR * C_LIGHT
V_HIGGS = 1.0 / np.sqrt(np.sqrt(2.0) * G_F)

M_TOP = 172.76
V_TB, V_TS = 0.99915, 0.0404
V_TD = 0.0086

M_B_QUARK = 4.18
M_S_QUARK = 0.095
M_C_QUARK = 1.30
M_D_QUARK = 0.0047

_M = lambda pid: Particle.from_pdgid(pid).mass * 1e-3
_TAU = lambda pid: Particle.from_pdgid(pid).lifetime * 1e-9

M_BPLUS = _M(521)
M_B0 = _M(511)
M_BS = _M(531)
M_LAMBDA_B = _M(5122)
M_KPLUS = _M(321)
M_K0 = _M(311)
M_LAMBDA = _M(3122)
M_PIPLUS = _M(211)
M_PI0 = _M(111)
M_D0 = _M(421)
M_ELECTRON = _M(11)
M_MUON = _M(13)
M_TAU = _M(15)

TAU_BPLUS = _TAU(521)
TAU_B0 = _TAU(511)
TAU_BS = _TAU(531)
TAU_LAMBDA_B = _TAU(5122)
TAU_KPLUS = _TAU(321)

M_S_MAX_BTOK = M_BPLUS - M_KPLUS

M_SPECTATOR = 2.0
C_4PI = 5.1e-9

ALPHA_S = 0.30


_WINKLER_CSV = Path(__file__).resolve().parent / "data" / "winkler_widths.csv"


def _load_winkler(path=_WINKLER_CSV):
    import csv
    grids = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            g = float(row["gamma_GeV"])
            if g > 0.0:
                grids.setdefault(row["curve"], []).append(
                    (float(row["m_phi_GeV"]), g))
    out = {}
    for curve, pts in grids.items():
        pts.sort()
        m = np.array([p[0] for p in pts])
        gm = np.array([p[1] for p in pts])
        out[curve] = (np.log(m), np.log(gm))
    return out


_WINKLER = _load_winkler()


def _winkler_width(curve, m_S):
    """Winkler Fig."""
    lm, lg = _WINKLER[curve]
    x = np.log(float(m_S))
    if x < lm[0] or x > lm[-1]:
        return 0.0
    return float(np.exp(np.interp(x, lm, lg)))


def _beta(m_S, m_d):
    """Two-body velocity sqrt(1 - 4 m_d^2 / m_S^2); 0 below threshold."""
    m_S = np.asarray(m_S, float)
    out = 1.0 - 4.0 * m_d ** 2 / m_S ** 2
    return np.sqrt(np.clip(out, 0.0, None))


def width_leptonic(m_S, m_lep):
    """Gamma(S -> l+ l-) at s_theta^2 = 1 (Winkler eq."""
    m_S = np.asarray(m_S, float)
    beta = _beta(m_S, m_lep)
    return G_F * m_S / (4.0 * np.sqrt(2.0) * np.pi) * m_lep ** 2 * beta ** 3


def _chpt_form_factors(m_S):
    """Leading-order ChPT scalar form factors (Winkler eq."""
    s = np.asarray(m_S, float) ** 2
    amp_pi = (7.0 / 9.0) * M_PIPLUS ** 2 + (2.0 / 9.0) * (s + 2.0 * M_PIPLUS ** 2)
    amp_K = ((7.0 / 9.0) * (0.5 * M_PIPLUS ** 2 + (M_KPLUS ** 2 - 0.5 * M_PIPLUS ** 2))
             + (2.0 / 9.0) * (s + 2.0 * M_KPLUS ** 2))
    return amp_pi, amp_K


def width_pipi(m_S):
    """Gamma(S -> pi pi) at s_theta^2 = 1, LO-ChPT form factors (eq."""
    m_S = np.asarray(m_S, float)
    beta = _beta(m_S, M_PIPLUS)
    amp_pi, _ = _chpt_form_factors(m_S)
    pref = 3.0 * G_F / (16.0 * np.sqrt(2.0) * np.pi * m_S)
    return np.where(beta > 0, pref * beta * amp_pi ** 2, 0.0)


def width_KK(m_S):
    """Gamma(S -> K K) at s_theta^2 = 1, LO-ChPT form factors (eq."""
    m_S = np.asarray(m_S, float)
    beta = _beta(m_S, M_KPLUS)
    _, amp_K = _chpt_form_factors(m_S)
    pref = G_F / (4.0 * np.sqrt(2.0) * np.pi * m_S)
    return np.where(beta > 0, pref * beta * amp_K ** 2, 0.0)


def width_4pi(m_S):
    """Gamma(S -> 4 pi, eta eta, rho rho, ...) at s_theta^2 = 1 (eq."""
    m_S = np.asarray(m_S, float)
    beta = _beta(m_S, 2.0 * M_PIPLUS)
    return np.where(beta > 0, C_4PI * m_S ** 3 * beta, 0.0)


def width_qq(m_S, m_quark, m_threshold_meson):
    """Spectator Gamma(S -> q qbar) at s_theta^2 = 1 (eq."""
    m_S = np.asarray(m_S, float)
    beta = _beta(m_S, m_threshold_meson)
    return G_F * m_S / (4.0 * np.sqrt(2.0) * np.pi) * 3.0 * m_quark ** 2 * beta ** 3


def width_gg(m_S):
    """Loop-induced Gamma(S -> g g) at s_theta^2 = 1 (eqs."""
    m_S = np.asarray(m_S, float)

    def f(x):
        x = np.asarray(x, float)
        out = np.empty(x.shape, dtype=complex)
        below = x <= 1.0
        out[below] = np.arcsin(np.sqrt(np.clip(x[below], 0, 1))) ** 2
        xa = x[~below]
        if xa.size:
            r = np.sqrt(1.0 - 1.0 / xa)
            out[~below] = -0.25 * (np.log((1.0 + r) / (1.0 - r)) - 1j * np.pi) ** 2
        return out

    amp = np.zeros(m_S.shape, dtype=complex)
    for m_q in (M_C_QUARK, M_B_QUARK, M_TOP):
        x = m_S ** 2 / (4.0 * m_q ** 2)
        amp = amp + (x + (x - 1.0) * f(x)) / x ** 2
    return (ALPHA_S ** 2 * m_S ** 3) / (32.0 * np.pi ** 3 * V_HIGGS ** 2) * np.abs(amp) ** 2


_LEPTONS = {"ee": M_ELECTRON, "mumu": M_MUON, "tautau": M_TAU}


WIDTH_SCHEMES = ("winkler", "chpt_spectator")


def partial_widths(m_S, scheme="winkler"):
    """All partial widths at ``s_theta^2 = 1`` for a scalar mass ``m_S`` (GeV)."""
    if scheme not in WIDTH_SCHEMES:
        raise ValueError(f"unknown width scheme {scheme!r}; expected one of {WIDTH_SCHEMES}")
    m_S = float(m_S)
    w = {}
    for name, m_l in _LEPTONS.items():
        w[name] = float(width_leptonic(m_S, m_l))
    if scheme == "chpt_spectator":
        if m_S < M_SPECTATOR:
            w["pipi"] = float(width_pipi(m_S))
            w["KK"] = float(width_KK(m_S))
            w["4pi"] = float(width_4pi(m_S))
            w["ss"] = w["cc"] = w["gg"] = 0.0
        else:
            w["pipi"] = w["KK"] = w["4pi"] = 0.0
            w["ss"] = float(width_qq(m_S, M_S_QUARK, M_KPLUS))
            w["cc"] = float(width_qq(m_S, M_C_QUARK, M_D0))
            w["gg"] = float(width_gg(m_S))
        return w

    if m_S < M_SPECTATOR:
        w["pipi"] = _winkler_width("pipi", m_S)
        w["KK"] = _winkler_width("KK", m_S)
        w["4pi"] = _winkler_width("4pi_etaeta_rhorho_extra", m_S)
        w["ss"] = w["cc"] = w["gg"] = 0.0
    else:
        w["pipi"] = w["KK"] = w["4pi"] = 0.0
        w["ss"] = _winkler_width("ss", m_S)
        w["cc"] = _winkler_width("cc", m_S)
        w["gg"] = _winkler_width("gg", m_S)
    return w


def total_width(m_S, scheme="winkler"):
    """Total width at ``s_theta^2 = 1`` (GeV)."""
    return sum(partial_widths(m_S, scheme=scheme).values())


def branching_ratios(m_S, scheme="winkler"):
    """Visible branching ratios at ``m_S`` (theta-independent: ratios of widths)."""
    w = partial_widths(m_S, scheme=scheme)
    tot = sum(w.values())
    if tot <= 0:
        return {k: 0.0 for k in w}
    return {k: v / tot for k, v in w.items()}


def ctau(m_S, sin2theta, scheme="winkler"):
    """Proper decay length c*tau in metres for mass ``m_S`` and coupling ``sin^2 theta``."""
    g = total_width(m_S, scheme=scheme) * float(sin2theta)
    if g <= 0:
        return np.inf
    return HBAR_C / g


def ctau_sin2theta1(m_S, scheme="winkler"):
    """c*tau (m) at ``sin^2 theta = 1`` -- the lifetime the scan divides by the coupling, the BC4
    analogue of HNLCalc's ``ctau(U^2 = 1)``."""
    g = total_width(m_S, scheme=scheme)
    return HBAR_C / g if g > 0 else np.inf


def g_phisb(sin2theta=1.0):
    """Effective b-s-S coupling g_{phi s b} (dimensionless, Winkler eq."""
    s_theta = np.sqrt(float(sin2theta))
    return (s_theta * M_B_QUARK / V_HIGGS) * (
        3.0 * np.sqrt(2.0) * G_F * M_TOP ** 2 * V_TS * V_TB) / (16.0 * np.pi ** 2)


def _lam(mx, my, mz):
    """Dimensionless Kallen factor of Winkler eq."""
    return ((mx ** 2 - (my - mz) ** 2) * (mx ** 2 - (my + mz) ** 2)) / mx ** 4


def _f_K(q2):
    """B -> K scalar form factor f_K(q^2) (Winkler eq."""
    return 0.33 / (1.0 - q2 / 37.5)


def width_B_to_K_S(m_S, m_B, sin2theta=1.0):
    """Gamma(B -> K S) (GeV) for parent mass ``m_B`` (Winkler eqs."""
    m_S = np.asarray(m_S, float)
    g2 = g_phisb(sin2theta) ** 2
    matrix2 = 0.25 * (m_B ** 2 - M_KPLUS ** 2) ** 2 / (M_B_QUARK - M_S_QUARK) ** 2 * _f_K(m_S ** 2) ** 2
    lam = _lam(m_B, M_KPLUS, m_S)
    lam_half = np.sqrt(np.clip(lam, 0.0, None))
    return g2 * matrix2 * lam_half / (16.0 * np.pi * m_B)


def br_B_to_K_S(m_S, parent="B+", sin2theta=1.0):
    """BR(B -> K S) for parent ``B+``/``B0``."""
    if parent == "B+":
        m_B, tau_B = M_BPLUS, TAU_BPLUS
    elif parent == "B0":
        m_B, tau_B = M_B0, TAU_B0
    else:
        raise ValueError(f"parent must be 'B+' or 'B0', got {parent!r}")
    if np.any(np.asarray(m_S) >= m_B - M_KPLUS):
        return np.where(np.asarray(m_S) >= m_B - M_KPLUS, 0.0,
                        width_B_to_K_S(m_S, m_B, sin2theta) * tau_B / HBAR)
    return width_B_to_K_S(m_S, m_B, sin2theta) * tau_B / HBAR


def br_B_to_Xs_S(m_S, parent="B+", sin2theta=1.0):
    """Inclusive BR(b-hadron -> X_s S), spectator estimate (Winkler eq."""
    m_B, tau_B = {"B+": (M_BPLUS, TAU_BPLUS),
                  "B0": (M_B0, TAU_B0),
                  "Bs": (M_BS, TAU_BS),
                  "Lambda_b": (M_LAMBDA_B, TAU_LAMBDA_B)}[parent]
    m_S = np.asarray(m_S, float)
    g2 = g_phisb(sin2theta) ** 2
    gamma = np.where(m_S < m_B,
                     g2 * (m_B ** 2 - m_S ** 2) ** 2 / (32.0 * np.pi * m_B ** 3), 0.0)
    return gamma * tau_B / HBAR


def g_phids(sin2theta=1.0):
    """Effective s-d-S coupling, eq."""
    s_theta = np.sqrt(float(sin2theta))
    return (s_theta * M_S_QUARK / V_HIGGS) * (
        3.0 * np.sqrt(2.0) * G_F * M_TOP ** 2 * V_TD * V_TS) / (16.0 * np.pi ** 2)


def br_K_to_pi_S(m_S, sin2theta=1.0):
    """BR(K+ -> pi+ S) (Winkler eqs."""
    m_S = np.asarray(m_S, float)
    g2 = g_phids(sin2theta) ** 2
    matrix2 = (0.5 * (M_KPLUS ** 2 - M_PIPLUS ** 2) / (M_S_QUARK - M_D_QUARK)) ** 2
    lam = _lam(M_KPLUS, M_PIPLUS, m_S)
    lam_half = np.sqrt(np.clip(lam, 0.0, None))
    gamma = np.where(m_S < M_KPLUS - M_PIPLUS,
                     g2 * matrix2 * lam_half / (16.0 * np.pi * M_KPLUS), 0.0)
    return gamma * TAU_KPLUS / HBAR
