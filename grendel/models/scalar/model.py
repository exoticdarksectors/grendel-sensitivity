"""Scalar-portal model layer (PBC benchmark BC4).

The light dark scalar ``S`` mixes with the Higgs through a single angle ``theta``
(``s_theta = sin theta``).  Below the electroweak scale it behaves as a Higgs
with universally suppressed couplings ``L = -(s_theta m_f / v) S f f``
(Winkler arXiv:1809.01876 eq. 2).  This module is the BC4 analogue of the HNL's
``HNLCalc``: it maps the single coupling ``sin^2 theta`` (at fixed ``m_S``) onto
the three quantities the sensitivity needs SIMULTANEOUSLY:

  * production rate   -> inclusive BR(B -> X_s S)(theta)
                         [with exclusive B -> K S and K -> pi S helpers]
  * lifetime          -> c*tau(m_S, theta)
  * visible decay     -> branching ratios S -> mu mu, ee, tau tau, pi pi, K K,
                         s s, c c, g g vs mass

Because BR(B->K S), every partial width, and hence c*tau all scale as
``s_theta^2``, while the *visible branching fractions are ratios of widths and so
are theta-independent, the coupling factors out exactly as ``U^2`` does for the
HNL.  Every width below is therefore returned at ``s_theta^2 = 1`` and the
sensitivity scan multiplies / divides by the actual ``sin^2 theta`` (mirroring
``grendel.reco.acceptance.scan_u2``).

References
----------
* M. W. Winkler, "Decay and Detection of a Light Scalar Boson Mixing with the
  Higgs", Phys. Rev. D 99 (2019) 015018, arXiv:1809.01876.  Equation numbers
  below (eq. 12, 15, 21, 30-33, A2-A9) refer to this paper.
* Unified FIP sensitivity calculation, arXiv:2311.00507 (the source the 2025 PBC
  report arXiv:2505.00947 uses for BC4 conventions and competitor curves): the
  BC4 parameter plane is ``(m_S, sin^2 theta)`` with B-meson production, the
  same Winkler decay rates, and ``N_signal >= 3`` background-free limits. The
  canonical GRENDEL normalization is the inclusive ``b -> X_s S`` rate; a
  two-body ``B -> K S`` recoil supplies the unobserved-system kinematics.

The central result uses digitized matched/dispersive widths from Winkler Fig. 4.
The analytic leading-order ChPT form factors below 2 GeV and perturbative
spectator widths above 2 GeV remain available as the ``chpt_spectator`` scheme;
their contour displacement defines the BC4 decay-model uncertainty envelope.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from particle import Particle

# ---------------------------------------------------------------------------
# Constants (PDG; same `particle` source the HNL production constants use)
# ---------------------------------------------------------------------------
G_F = 1.1663787e-5                     # Fermi constant, GeV^-2
HBAR = 6.582119569e-25                 # reduced Planck constant, GeV*s
C_LIGHT = 299792458.0                  # m/s
HBAR_C = HBAR * C_LIGHT                # 1.973269804e-16 GeV*m
V_HIGGS = 1.0 / np.sqrt(np.sqrt(2.0) * G_F)   # Higgs vev, ~246.22 GeV

# Loop / CKM inputs entering the b -> s S (and s -> d S) effective coupling.
M_TOP = 172.76                         # GeV (PDG pole mass)
V_TB, V_TS = 0.99915, 0.0404           # |V_tb|, |V_ts|
V_TD = 0.0086                          # |V_td|  (for K -> pi S)

# Quark masses entering the effective couplings / spectator widths.
# MSbar for b; Winkler footnote 9 uses m_s = 95 MeV, m_c = 1.3 GeV.
M_B_QUARK = 4.18
M_S_QUARK = 0.095
M_C_QUARK = 1.30
M_D_QUARK = 0.0047

_M = lambda pid: Particle.from_pdgid(pid).mass * 1e-3       # GeV
_TAU = lambda pid: Particle.from_pdgid(pid).lifetime * 1e-9  # s

M_BPLUS = _M(521)
M_B0 = _M(511)
M_BS = _M(531)
M_LAMBDA_B = _M(5122)   # lightest b-baryon; proxy for the whole b-baryon pool
M_KPLUS = _M(321)
M_K0 = _M(311)
M_LAMBDA = _M(3122)     # lightest s-baryon: the Lambda_b -> Lambda S recoil,
                        # the baryonic analogue of the kaon recoil for mesons
M_PIPLUS = _M(211)
M_PI0 = _M(111)
M_D0 = _M(421)        # lightest c-meson: spectator c-cbar threshold
M_ELECTRON = _M(11)
M_MUON = _M(13)
M_TAU = _M(15)

TAU_BPLUS = _TAU(521)
TAU_B0 = _TAU(511)
TAU_BS = _TAU(531)
TAU_LAMBDA_B = _TAU(5122)
TAU_KPLUS = _TAU(321)

# B -> K S kinematic ceiling (production closes here).
M_S_MAX_BTOK = M_BPLUS - M_KPLUS       # ~4.785 GeV

# Hadronic-region matching point: dispersive (pi pi / K K) below, perturbative
# spectator (s s, c c, g g) above (Winkler sec. III D, valid for m_S > 2 GeV).
M_SPECTATOR = 2.0                      # GeV
# 4 pi, eta eta, rho rho contribution (Winkler eq. 33), tuned so the hadronic
# rate transits smoothly into the spectator model at 2 GeV.
C_4PI = 5.1e-9                         # GeV^-2

# alpha_s(m_S) for the loop-induced g g width (eq. 31).  A slow function over
# the 0.5-5 GeV window; the central value below is adequate at the ~10% level
# the spectator model itself carries.
ALPHA_S = 0.30


# ---------------------------------------------------------------------------
# Winkler dispersive hadronic widths (digitized from arXiv:1809.01876 Fig. 4;
# see scalar/data/winkler_widths.csv + its extractor README). These supersede
# the LO-ChPT sector below 2 GeV and the perturbative-spectator widths above
# with the paper's *matched* result -- fixing the missing f0(980) enhancement
# and the ~factor-few seam at 2 GeV. Per-channel curves are kept separate so
# the acceptance's decay-mode breakdown still works. All at s_theta^2 = 1;
# interpolation is linear in log(m)-log(Gamma).
# ---------------------------------------------------------------------------
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
    """Winkler Fig. 4 width Gamma_curve(m_S) at s_theta^2 = 1 (GeV); 0 outside
    the curve's digitized mass range."""
    lm, lg = _WINKLER[curve]
    x = np.log(float(m_S))
    if x < lm[0] or x > lm[-1]:
        return 0.0
    return float(np.exp(np.interp(x, lm, lg)))


# ---------------------------------------------------------------------------
# 1. Visible decay widths  (all at s_theta^2 = 1)
# ---------------------------------------------------------------------------
def _beta(m_S, m_d):
    """Two-body velocity sqrt(1 - 4 m_d^2 / m_S^2); 0 below threshold."""
    m_S = np.asarray(m_S, float)
    out = 1.0 - 4.0 * m_d ** 2 / m_S ** 2
    return np.sqrt(np.clip(out, 0.0, None))


def width_leptonic(m_S, m_lep):
    """Gamma(S -> l+ l-) at s_theta^2 = 1  (Winkler eq. 12)."""
    m_S = np.asarray(m_S, float)
    beta = _beta(m_S, m_lep)
    return G_F * m_S / (4.0 * np.sqrt(2.0) * np.pi) * m_lep ** 2 * beta ** 3


def _chpt_form_factors(m_S):
    """Leading-order ChPT scalar form factors (Winkler eq. 21) at sqrt(s)=m_S,
    contracted with the scalar coupling (7/9)(Gamma + Delta) + (2/9) Theta of
    his eq. 13.  Returns the pi pi and K K amplitudes (GeV^2)."""
    s = np.asarray(m_S, float) ** 2
    # pi pi : Gamma_pi = m_pi^2, Delta_pi = 0, Theta_pi = s + 2 m_pi^2
    amp_pi = (7.0 / 9.0) * M_PIPLUS ** 2 + (2.0 / 9.0) * (s + 2.0 * M_PIPLUS ** 2)
    # K K : Gamma_K = m_pi^2/2, Delta_K = m_K^2 - m_pi^2/2, Theta_K = s + 2 m_K^2
    amp_K = ((7.0 / 9.0) * (0.5 * M_PIPLUS ** 2 + (M_KPLUS ** 2 - 0.5 * M_PIPLUS ** 2))
             + (2.0 / 9.0) * (s + 2.0 * M_KPLUS ** 2))
    return amp_pi, amp_K


def width_pipi(m_S):
    """Gamma(S -> pi pi) at s_theta^2 = 1, LO-ChPT form factors (eq. 15 + 21)."""
    m_S = np.asarray(m_S, float)
    beta = _beta(m_S, M_PIPLUS)
    amp_pi, _ = _chpt_form_factors(m_S)
    pref = 3.0 * G_F / (16.0 * np.sqrt(2.0) * np.pi * m_S)
    return np.where(beta > 0, pref * beta * amp_pi ** 2, 0.0)


def width_KK(m_S):
    """Gamma(S -> K K) at s_theta^2 = 1, LO-ChPT form factors (eq. 15 + 21)."""
    m_S = np.asarray(m_S, float)
    beta = _beta(m_S, M_KPLUS)
    _, amp_K = _chpt_form_factors(m_S)
    pref = G_F / (4.0 * np.sqrt(2.0) * np.pi * m_S)
    return np.where(beta > 0, pref * beta * amp_K ** 2, 0.0)


def width_4pi(m_S):
    """Gamma(S -> 4 pi, eta eta, rho rho, ...) at s_theta^2 = 1 (eq. 33)."""
    m_S = np.asarray(m_S, float)
    beta = _beta(m_S, 2.0 * M_PIPLUS)   # 4 pi threshold via beta_{2pi}
    return np.where(beta > 0, C_4PI * m_S ** 3 * beta, 0.0)


def width_qq(m_S, m_quark, m_threshold_meson):
    """Spectator Gamma(S -> q qbar) at s_theta^2 = 1 (eq. 30): the leptonic
    prefactor with the colour factor 3 and the quark mass, the velocity set by
    the lightest meson carrying that quark (Winkler footnote 9)."""
    m_S = np.asarray(m_S, float)
    beta = _beta(m_S, m_threshold_meson)
    return G_F * m_S / (4.0 * np.sqrt(2.0) * np.pi) * 3.0 * m_quark ** 2 * beta ** 3


def width_gg(m_S):
    """Loop-induced Gamma(S -> g g) at s_theta^2 = 1 (eqs. 31-32)."""
    m_S = np.asarray(m_S, float)

    def f(x):
        # x > 1 (m_S > 2 m_q) gives a complex f(x); |amp|^2 is taken below, so
        # the array must stay complex to retain the -i pi piece (Winkler eq. 32).
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


# Leptonic / partonic final states tracked for the visible-BR table.
_LEPTONS = {"ee": M_ELECTRON, "mumu": M_MUON, "tautau": M_TAU}


WIDTH_SCHEMES = ("winkler", "chpt_spectator")


def partial_widths(m_S, scheme="winkler"):
    """All partial widths at ``s_theta^2 = 1`` for a scalar mass ``m_S`` (GeV).

    Returns an ordered dict {channel: Gamma_GeV}.  Below ``M_SPECTATOR`` the
    hadronic sector is dispersive (pi pi, K K, 4 pi); above it the perturbative
    spectator model (s s, c c, g g).  The cross-over is the Winkler prescription
    (his Fig. 4); it is not perfectly continuous at 2 GeV (the eq.-33 constant is
    tuned to soften the seam) but is smooth at the ~order-unity level.
    """
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

    # Hadronic sector from Winkler's matched dispersive result (Fig. 4),
    # switching at M_SPECTATOR = 2 GeV (his matching point) between the
    # dispersive pi pi / K K / 4 pi channels below and the perturbative
    # s s / c c / g g continuum above. The digitized per-channel curves meet
    # continuously at 2 GeV; gating avoids the narrow double-count where the
    # two sets overlap at the boundary. Supersedes the analytic width_pipi/KK/
    # 4pi (LO ChPT) and width_qq/gg (spectator), kept above as cross-check.
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
    """Total width at ``s_theta^2 = 1`` (GeV).  Gamma_tot(theta) = s_theta^2 * this."""
    return sum(partial_widths(m_S, scheme=scheme).values())


def branching_ratios(m_S, scheme="winkler"):
    """Visible branching ratios at ``m_S`` (theta-independent: ratios of widths)."""
    w = partial_widths(m_S, scheme=scheme)
    tot = sum(w.values())
    if tot <= 0:
        return {k: 0.0 for k in w}
    return {k: v / tot for k, v in w.items()}


def ctau(m_S, sin2theta, scheme="winkler"):
    """Proper decay length c*tau in metres for mass ``m_S`` and coupling
    ``sin^2 theta``.  c*tau = hbar c / (sin^2 theta * Gamma_hat(m_S))."""
    g = total_width(m_S, scheme=scheme) * float(sin2theta)
    if g <= 0:
        return np.inf
    return HBAR_C / g


def ctau_sin2theta1(m_S, scheme="winkler"):
    """c*tau (m) at ``sin^2 theta = 1`` -- the lifetime the scan divides by the
    coupling, the BC4 analogue of HNLCalc's ``ctau(U^2 = 1)``."""
    g = total_width(m_S, scheme=scheme)
    return HBAR_C / g if g > 0 else np.inf


# ---------------------------------------------------------------------------
# 2. Production:  B -> K S  (and inclusive B -> X_s S, and K -> pi S)
# ---------------------------------------------------------------------------
def g_phisb(sin2theta=1.0):
    """Effective b-s-S coupling g_{phi s b} (dimensionless, Winkler eq. A2),
    scaled to ``sin^2 theta`` (it is linear in s_theta, so |g|^2 ~ sin^2 theta)."""
    s_theta = np.sqrt(float(sin2theta))
    return (s_theta * M_B_QUARK / V_HIGGS) * (
        3.0 * np.sqrt(2.0) * G_F * M_TOP ** 2 * V_TS * V_TB) / (16.0 * np.pi ** 2)


def _lam(mx, my, mz):
    """Dimensionless Kallen factor of Winkler eq. A4."""
    return ((mx ** 2 - (my - mz) ** 2) * (mx ** 2 - (my + mz) ** 2)) / mx ** 4


def _f_K(q2):
    """B -> K scalar form factor f_K(q^2) (Winkler eq. A6), q^2 = m_S^2."""
    return 0.33 / (1.0 - q2 / 37.5)


def width_B_to_K_S(m_S, m_B, sin2theta=1.0):
    """Gamma(B -> K S) (GeV) for parent mass ``m_B`` (Winkler eqs. A3, A5, A6)."""
    m_S = np.asarray(m_S, float)
    g2 = g_phisb(sin2theta) ** 2
    matrix2 = 0.25 * (m_B ** 2 - M_KPLUS ** 2) ** 2 / (M_B_QUARK - M_S_QUARK) ** 2 * _f_K(m_S ** 2) ** 2
    lam = _lam(m_B, M_KPLUS, m_S)
    lam_half = np.sqrt(np.clip(lam, 0.0, None))
    return g2 * matrix2 * lam_half / (16.0 * np.pi * m_B)


def br_B_to_K_S(m_S, parent="B+", sin2theta=1.0):
    """BR(B -> K S) for parent ``B+``/``B0``.  Linear in ``sin^2 theta``; the
    production driver evaluates this at sin2theta=1 and the scan multiplies."""
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
    """Inclusive BR(b-hadron -> X_s S), spectator estimate (Winkler eq. A7);
    ~10x the exclusive K mode at low mass (~5.3 sin^2 theta).  This is the
    production normalization used by the scan -- GRENDEL reconstructs only the
    S vertex, so the prompt X_s system is summed over inclusively.  g_{phi s b}
    is a b-quark process (spectator-independent), so only m_B, tau vary -- which
    is exactly why b-baryons enter here on the same footing as the mesons, with
    just the Lambda_b mass and lifetime (the b-baryon pool proxy)."""
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
    """Effective s-d-S coupling, eq. A2 under (b,s)->(s,d)."""
    s_theta = np.sqrt(float(sin2theta))
    return (s_theta * M_S_QUARK / V_HIGGS) * (
        3.0 * np.sqrt(2.0) * G_F * M_TOP ** 2 * V_TD * V_TS) / (16.0 * np.pi ** 2)


def br_K_to_pi_S(m_S, sin2theta=1.0):
    """BR(K+ -> pi+ S) (Winkler eqs. A8, A9); relevant for m_S < m_K - m_pi.

    Implemented for model completeness (the very-low-mass production channel).
    Sampling a kaon flux at the LHC IP is not implemented.
    """
    m_S = np.asarray(m_S, float)
    g2 = g_phids(sin2theta) ** 2
    matrix2 = (0.5 * (M_KPLUS ** 2 - M_PIPLUS ** 2) / (M_S_QUARK - M_D_QUARK)) ** 2
    lam = _lam(M_KPLUS, M_PIPLUS, m_S)
    lam_half = np.sqrt(np.clip(lam, 0.0, None))
    gamma = np.where(m_S < M_KPLUS - M_PIPLUS,
                     g2 * matrix2 * lam_half / (16.0 * np.pi * M_KPLUS), 0.0)
    return gamma * TAU_KPLUS / HBAR
