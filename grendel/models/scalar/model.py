"""Scalar-portal model layer (PBC benchmarks BC4 and BC5)."""
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

ALPHA_S_MZ = 0.1180
M_Z = 91.1876
_ZETA3 = 1.2020569031595942


def _qcd_beta(a, nf):
    """-da/dln(mu^2) for a = alpha_s/pi, four loops (MSbar)."""
    b0 = (11.0 - 2.0 / 3.0 * nf) / 4.0
    b1 = (102.0 - 38.0 / 3.0 * nf) / 16.0
    b2 = (2857.0 / 2.0 - 5033.0 / 18.0 * nf + 325.0 / 54.0 * nf ** 2) / 64.0
    b3 = ((149753.0 / 6.0 + 3564.0 * _ZETA3) - (1078361.0 / 162.0 + 6508.0 / 27.0 * _ZETA3) * nf
          + (50065.0 / 162.0 + 6472.0 / 81.0 * _ZETA3) * nf ** 2 + 1093.0 / 729.0 * nf ** 3) / 256.0
    return a * a * (b0 + a * (b1 + a * (b2 + a * b3)))


def _run(a, mu_from, mu_to, nf, steps=200):
    """RK4 in ln(mu^2) from mu_from to mu_to at fixed nf."""
    h = 2.0 * np.log(mu_to / mu_from) / steps
    for _ in range(steps):
        k1 = -_qcd_beta(a, nf)
        k2 = -_qcd_beta(a + 0.5 * h * k1, nf)
        k3 = -_qcd_beta(a + 0.5 * h * k2, nf)
        k4 = -_qcd_beta(a + h * k3, nf)
        a = a + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return a


def alpha_s(mu):
    """alpha_s(mu) for scalar or array ``mu`` (GeV), run down from M_Z."""
    mu = np.asarray(mu, float)
    out = np.empty(mu.shape)
    for i, m in np.ndenumerate(mu):
        a = _run(ALPHA_S_MZ / np.pi, M_Z, max(m, M_B_QUARK), 5)
        if m < M_B_QUARK:
            a = _run(a, M_B_QUARK, max(m, M_C_QUARK), 4)
        if m < M_C_QUARK:
            a = _run(a, M_C_QUARK, m, 3)
        out[i] = np.pi * a
    return out if out.ndim else float(out)


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

WINKLER_M_MAX = 7.5
_EVENTCALC_DIR = Path(__file__).resolve().parent / "data" / "eventcalc_1809"
_EVENTCALC_CHANNELS = {"ePeM": "ee", "muPmuM": "mumu", "tauPtauM": "tautau",
                       "Jets-ss": "ss", "Jets-cc": "cc", "Jets-GG": "gg",
                       "Jets-bb": "bb"}


def _load_eventcalc(directory=_EVENTCALC_DIR):
    import gzip
    import json
    with gzip.open(directory / "ctau.json.gz", "rt") as fh:
        ctau = np.array(json.load(fh), float)
    with gzip.open(directory / "branching_ratios.json.gz", "rt") as fh:
        channels = json.load(fh)
    br = {}
    for entry in channels:
        pts = np.array(entry[2], float)
        br[entry[0]] = (pts[:, 0], pts[:, 1])
    return (np.log(ctau[:, 0]), np.log(ctau[:, 1])), br


_EVENTCALC = None


def _eventcalc():
    global _EVENTCALC
    if _EVENTCALC is None:
        _EVENTCALC = _load_eventcalc()
    return _EVENTCALC


def eventcalc_ctau_sin2theta1(m_S):
    """EventCalc c*tau (m) at sin^2 theta = 1, log-log interpolated (0.01-63 GeV; inf outside)."""
    lm, lc = _eventcalc()[0]
    x = np.log(float(m_S))
    if x < lm[0] or x > lm[-1]:
        return np.inf
    return float(np.exp(np.interp(x, lm, lc)))


def eventcalc_branching_ratios(m_S):
    """EventCalc branching ratios {channel: BR} at ``m_S`` (their 16 channels)."""
    m_S = float(m_S)
    return {name: float(np.interp(m_S, m, b, left=0.0, right=0.0))
            for name, (m, b) in _eventcalc()[1].items()}


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
    return (alpha_s(m_S) ** 2 * m_S ** 3) / (32.0 * np.pi ** 3 * V_HIGGS ** 2) * np.abs(amp) ** 2


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
        w["bb"] = float(width_qq(m_S, M_B_QUARK, M_BPLUS)) if m_S >= M_SPECTATOR else 0.0
        return w

    if m_S > WINKLER_M_MAX:
        total = HBAR_C / eventcalc_ctau_sin2theta1(m_S)
        br = eventcalc_branching_ratios(m_S)
        mapped = sum(br[name] for name in _EVENTCALC_CHANNELS)
        if not np.isfinite(total) or mapped < 0.98:
            raise ValueError(f"EventCalc tables do not cover m_S = {m_S} GeV "
                             f"(mapped BR {mapped:.3f})")
        w = {k: 0.0 for k in ("ee", "mumu", "tautau", "pipi", "KK", "4pi", "ss", "cc", "gg", "bb")}
        for name, key in _EVENTCALC_CHANNELS.items():
            w[key] = total * br[name] / mapped
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
    w["bb"] = 0.0
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
    """BR(B -> K S) for parent ``B+``/``B0``: the exclusive reference; production uses
    ``br_B_to_Xs_S`` (~11x larger at 0.5 GeV)."""
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


M_HIGGS = 125.20
GAMMA_HIGGS_SM = 4.07e-3
XI_SB = 3.6e-4
BR_HSS_BC5 = 0.01
M_S_MAX_HSS = M_HIGGS / 2.0

_SENSCALC_DIR = Path(__file__).resolve().parent / "data" / "senscalc_quartic"
_SENSCALC_TABLES: dict = {}


def br_h_to_SS(m_S, alpha):
    """BR(h -> S S) for the quartic coupling ``alpha`` (GeV), Boiarska eq."""
    m_S = np.asarray(m_S, float)
    p_S = 0.5 * M_HIGGS * _beta(M_HIGGS, m_S)
    return alpha ** 2 * p_S / (16.0 * np.pi * M_HIGGS ** 2 * GAMMA_HIGGS_SM)


def alpha_quartic(m_S, br_hss=BR_HSS_BC5):
    """The quartic coupling ``alpha`` (GeV) that gives BR(h -> SS) = ``br_hss`` at ``m_S`` (inverse of
    :func:`br_h_to_SS`); 0 above m_h / 2."""
    p_S = 0.5 * M_HIGGS * float(_beta(M_HIGGS, float(m_S)))
    if p_S <= 0.0 or br_hss <= 0.0:
        return 0.0
    return float(np.sqrt(br_hss * 16.0 * np.pi * M_HIGGS ** 2 * GAMMA_HIGGS_SM / p_S))


def dgamma_dq2_B_to_K_SS(q2, m_S, alpha=1.0, m_B=M_BPLUS, m_recoil=M_KPLUS):
    """d Gamma / d q^2 (GeV^-1) of B -> K S S through the off-shell Higgs, Boiarska eq."""
    q2 = np.asarray(q2, float)
    q = np.sqrt(np.clip(q2, 0.0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        E2 = 0.5 * q
        E3 = (m_B ** 2 - q2 - m_recoil ** 2) / (2.0 * q)
        matrix = _f_K(q2) * (m_B ** 2 - M_KPLUS ** 2) / (M_B_QUARK - M_S_QUARK)
        phase = np.sqrt(np.clip(E2 ** 2 - m_S ** 2, 0.0, None)) * \
            np.sqrt(np.clip(E3 ** 2 - m_recoil ** 2, 0.0, None))
    pref = XI_SB ** 2 * M_B_QUARK ** 2 * alpha ** 2 / (512.0 * np.pi ** 3 * m_B ** 3 * V_HIGGS ** 2 * M_HIGGS ** 4)
    out = pref * matrix ** 2 * phase
    inside = (q2 > 4.0 * m_S ** 2) & (q2 < (m_B - m_recoil) ** 2)
    return np.where(inside, np.nan_to_num(out), 0.0)


def _q2_grid(m_S, m_B, m_recoil, n=400):
    lo, hi = 4.0 * m_S ** 2, (m_B - m_recoil) ** 2
    return np.linspace(lo, hi, n) if hi > lo else np.empty(0)


def br_B_to_K_SS(m_S, alpha=1.0, parent="B+"):
    """BR(B -> K S S) at the quartic coupling ``alpha`` (GeV): eq."""
    m_B, tau_B, m_K = {"B+": (M_BPLUS, TAU_BPLUS, M_KPLUS),
                       "B0": (M_B0, TAU_B0, M_K0)}[parent]
    q2 = _q2_grid(float(m_S), m_B, m_K)
    if len(q2) == 0:
        return 0.0
    gamma = float(np.trapezoid(dgamma_dq2_B_to_K_SS(q2, float(m_S), alpha, m_B, m_K), q2))
    return gamma * tau_B / HBAR


def sample_q2_B_to_K_SS(n, m_S, rng, m_B=M_BPLUS, m_recoil=M_KPLUS):
    """``n`` values of the S S invariant mass squared of B -> X_s S S drawn from eq."""
    q2 = _q2_grid(float(m_S), m_B, m_recoil)
    if len(q2) == 0:
        return np.empty(0)
    density = dgamma_dq2_B_to_K_SS(q2, float(m_S), 1.0, m_B, m_recoil)
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(q2))])
    if cdf[-1] <= 0.0:
        return np.empty(0)
    return np.interp(rng.random(n) * cdf[-1], cdf, q2)


def _senscalc_table(name):
    if name not in _SENSCALC_TABLES:
        data = np.loadtxt(_SENSCALC_DIR / f"{name}.dat")
        _SENSCALC_TABLES[name] = (data[:, 0], data[:, 1])
    return _SENSCALC_TABLES[name]


def _senscalc_br(name, m_S, alpha):
    m, br = _senscalc_table(name)
    m_S = np.asarray(m_S, float)
    return np.interp(m_S, m, br, left=br[0], right=0.0) * alpha ** 2


def br_b_to_Xs_SS(m_S, parent="B+", alpha=1.0):
    """Inclusive BR(b-hadron -> X_s S S) at the quartic coupling ``alpha`` (GeV): the SensCalc ``B+ ->
    X_s S S`` total, scaled by alpha^2 and, for the other b-hadrons, by the lifetime ratio -- b ->
    s S S is a b-quark process, so the parent enters only through its lifetime (and its recoil
    kinematics in the production driver), exactly as ``br_B_to_Xs_S``."""
    tau = {"B+": TAU_BPLUS, "B0": TAU_B0, "Bs": TAU_BS, "Lambda_b": TAU_LAMBDA_B}[parent]
    return _senscalc_br("BplustoXsSStotal", m_S, alpha) * tau / TAU_BPLUS


def br_Bs_to_SS(m_S, alpha=1.0):
    """BR(B_s -> S S) at the quartic coupling ``alpha`` (GeV), SensCalc table (Boiarska eq."""
    return _senscalc_br("BstoSS", m_S, alpha)
