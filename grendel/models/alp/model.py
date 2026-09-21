r"""BC10 fermiophilic-ALP model layer (couplings -> production, lifetime, BRs).

PBC benchmark **BC10**: a pseudoscalar axion-like particle ``a`` whose dominant
coupling is to SM fermions (the genuine charged-track "ALP" benchmark).  The
single parameter that this analysis sets a limit on is the inverse decay
constant ``1/f`` (universal Wilson coefficients ``c_f = 1``); it controls
production, lifetime, and the visible branching ratios *simultaneously*.

Coupling convention (fixes the meaning of the 1/f axis): the effective
pseudoscalar Yukawa g_aff = c_f m_f / f of Bauer-Neubert-Thamm
(arXiv:1708.00443), i.e.

    Gamma(a -> f fbar) = N_c (c_f m_f)^2 m_a (1/f)^2 beta / (8 pi).

The state-of-the-art BC10 phenomenology paper, Garcia-Kahlhoefer-Ovchynnikov-
Zaporozhchenko (GKOZ) arXiv:2310.03524, normalises its Lagrangian with
d_mu a/(2 f_GKOZ) instead, so our axis maps onto theirs as
1/f (here) = 2/f_GKOZ.  All external inputs below have been converted into the
BNT convention and cross-checked on physical (convention-free) anchors; see
data/senscalc_2501/PROVENANCE.md and data/alpinist/PROVENANCE.md.

External inputs (replacing the analytic placeholders of the first version):

1. **Decay widths and exclusive branching ratios (data-driven).**  The tables
   accompanying arXiv:2501.04525 are decoded from the exact SensCalc v.1.3.3
   MX inputs and committed as reviewable CSV/JSON data.  They cover
   0.01--10 GeV and include the eta/eta' resonance structure, the charm
   transition, and the exclusive low-energy hadronic modes.  The upstream
   convention Gamma=(gY/(2 vH))^2 coefficient maps to the BNT convention
   1/f=gY/vH by dividing every coefficient by four.  Direct anchors from the
   decoded table are BR(a->mumu)=0.202243 at exactly 1 GeV and the sharp
   charm-region structure around 2.58 GeV.

   The same SensCalc release also ships the exact arXiv:2310.03524 widths and
   branching functions; they define a one-sided heavy-pseudoscalar structural
   comparison, not a central prediction, and are not used here.

2. **b -> s a production coupling (one-loop RG, finite terms).**  The
   flavour-violating coefficient is evaluated with the ALPINIST implementation
   of the GKOZ RG equations (tools/compute_cbs_alpinist.py, ported under
   BSD-3):

       g_bs = CBS_EFF * (1/f),   CBS_EFF = 3.518383e-4  (Lambda_UV = 1 TeV)

   replacing the leading-log-with-C_TOP=1 placeholder (which was 2.5x larger
   in amplitude -- the finite one-loop and RG terms partially cancel the small
   log at 1 TeV).  Cross-check: GKOZ Table 1 quotes |C_bs| = 1.8e-3 at
   f_GKOZ = 1 GeV, i.e. c_bs = 2 |C_bs| / m_b ~ 7.5e-4 -- 8.6e-4 (m_b pole vs
   MSbar), whose BNT-convention half is 3.7e-4 -- 4.3e-4, within 6--20% of
   CBS_EFF.  The island's production-limited (lower) edge inherits this
   residual as <~20% in 1/f.

3. **Production channels: the full B -> K^(i) a tower** (K, K0*(700),
   K0*(1430), K*(892), K*(1410), K*(1680), K1(1270), K1(1400), K2*(1430)),
   for both B+ and B0, with the Boiarska et al. (arXiv:1904.10447) form
   factors as implemented in ALPINIST.  GKOZ find the tower gives ~4x the
   K + K*(892) rate alone.  (Bs -> phi a is not included, matching ALPINIST.)

The visible-channel classification follows SensCalc's ``procListnoecal``
selection from ``codes/EventCalc/DecayProductsSampler.nb``: gamma gamma,
3 pi0, and 2 K_L pi0 are excluded, while the remaining exclusive channels are
track-capable. This legacy classification remains available for comparison;
publication templates sample all exclusive modes in
generate_decay_templates_pythia.py and let reconstruction determine whether
the fully decayed stable final state contains reconstructable tracks.
"""

from __future__ import annotations

import json

import numpy as np
from particle import Particle

from ...io.tables import open_table, read_table_text
from .decay_models import DEFAULT_DECAY_MODEL, decay_data_dir

# --------------------------------------------------------------------------
# Physical constants (PDG 2024 central values)
# --------------------------------------------------------------------------
HBAR_C_GEV_M = 1.973269804e-16   # GeV * m  (hbar c)
HBAR_GEV_S = 6.582119569e-25     # GeV * s

# Lepton + quark masses (GeV).  Quark masses enter B-decay production below.
M_E = Particle.from_pdgid(11).mass * 1e-3
M_MU = Particle.from_pdgid(13).mass * 1e-3
M_TAU = Particle.from_pdgid(15).mass * 1e-3
M_PION = Particle.from_pdgid(211).mass * 1e-3

QUARK_MASSES = {                 # MS-bar current masses (GeV)
    "u": 0.00216, "d": 0.00467, "s": 0.0934, "c": 1.27, "b": 4.18,
}

# B masses and lifetimes for the B -> K^(i) a production rates.
M_BPLUS = Particle.from_pdgid(521).mass * 1e-3
M_B0 = Particle.from_pdgid(511).mass * 1e-3
TAU_BPLUS_S = 1.638e-12
TAU_B0_S = 1.519e-12

# --------------------------------------------------------------------------
# Production-coupling anchor (provenance: module docstring + data/alpinist/).
# Central decay inputs are loaded from data/senscalc_2501/ below.
# --------------------------------------------------------------------------
# b -> s a coefficient per unit (1/f): ALPINIST above_EW.C_qq_ij_mu at
# Lambda = 1 TeV, universal c_f = 1 (tools/compute_cbs_alpinist.py).
CBS_EFF = 3.518383e-4
# Reference inverse decay constant at which production weights and ctau are
# tabulated; the (m_a, 1/f) scan rescales off this point (production ~ (1/f)^2,
# ctau ~ (1/f)^-2), so the choice is just a numerical anchor.
INV_F_REF = 1.0e-3               # GeV^-1   (f_ref = 1 TeV)

# Light-meson poles where the perturbative ALP-meson diagonalization is not
# valid.  These are the exact ``exclRes`` windows used by the pinned SensCalc
# ALP-fermion analysis and quoted in arXiv:2501.04525 as excluded from its
# phenomenological interpretation.
LIGHT_MESON_RESONANCE_WINDOWS = (
    (0.125, 0.140, "pi0"),
    (0.538, 0.555, "eta"),
    (0.940, 0.974, "eta-prime"),
)

SENSCALC_2501_DATA_DIR = decay_data_dir(DEFAULT_DECAY_MODEL)

_WIDTH_CACHE: dict[str, tuple[np.ndarray, dict[str, np.ndarray]]] = {}
_BRANCHING_CACHE: dict[str, tuple[np.ndarray, dict[str, np.ndarray]]] = {}

# SensCalc's no-ECAL selection removes only these exclusive channels.  Channel
# IDs are stable products of the pinned exporter and are checked on load.
SENSCALC_INVISIBLE_CHANNEL_IDS = {
    "channel_004",  # gamma gamma
    "channel_014",  # 3 pi0
    "channel_021",  # 2 K_L pi0
}
# The MX list repeats the neutral and charged K* K* processes at positions
# 12/24 and 11/29 with byte-identical products and branching expressions.
# SensCalc installs them as DownValues keyed by process name, so the later
# definitions overwrite rather than add. Exclude the repeated positions to
# reproduce that behavior in the positional CSV representation.
SENSCALC_DUPLICATE_CHANNEL_IDS = {
    "channel_024",  # duplicate of channel_012
    "channel_029",  # duplicate of channel_011
}
SENSCALC_VISIBLE_HADRONIC_CHANNEL_IDS = tuple(
    f"channel_{index:03d}"
    for index in range(5, 33)
    if f"channel_{index:03d}" not in (
        SENSCALC_INVISIBLE_CHANNEL_IDS | SENSCALC_DUPLICATE_CHANNEL_IDS
    )
)


def excluded_light_meson_resonance(m_a):
    """Name the unsupported light-meson pole containing ``m_a``, if any."""
    mass = float(m_a)
    for lower, upper, name in LIGHT_MESON_RESONANCE_WINDOWS:
        if lower < mass < upper:
            return name
    return None


def _load_width_tables(decay_model=DEFAULT_DECAY_MODEL):
    """Return one pinned mass grid and its canonical BNT width columns."""
    if decay_model not in _WIDTH_CACHE:
        data_dir = decay_data_dir(decay_model)
        metadata = json.loads(
            (data_dir / "widths_metadata.json").read_text()
        )
        table = np.loadtxt(
            open_table(data_dir / "widths_bnt.csv"),
            delimiter=",",
            skiprows=1,
        )
        columns = {
            entry["canonical_name"]: table[:, entry["source_index"] - 1]
            for entry in metadata["columns"]
            if entry["canonical_name"] is not None
        }
        required = {
            "ee", "mumu", "tautau", "gammagamma",
            "nonhadronic_total", "hadronic_total", "total",
        }
        if not required <= columns.keys():
            missing = ", ".join(sorted(required - columns.keys()))
            raise RuntimeError(
                f"SensCalc {decay_model} width columns missing: {missing}"
            )
        _WIDTH_CACHE[decay_model] = (table[:, 0], columns)
    return _WIDTH_CACHE[decay_model]


def _load_branching_tables(decay_model=DEFAULT_DECAY_MODEL):
    """Return one pinned mass grid and exclusive branching-ratio columns."""
    if decay_model not in _BRANCHING_CACHE:
        data_dir = decay_data_dir(decay_model)
        channels = json.loads(
            read_table_text(data_dir / "decay_channels.json")
        )
        table = np.loadtxt(
            open_table(data_dir / "branching_ratios.csv"),
            delimiter=",",
            skiprows=1,
        )
        columns = {
            channel["id"]: table[:, index]
            for index, channel in enumerate(channels, start=1)
        }
        expected = {f"channel_{index:03d}" for index in range(1, 33)}
        if columns.keys() != expected:
            raise RuntimeError(
                f"SensCalc {decay_model} decay-channel IDs are incomplete"
            )
        _BRANCHING_CACHE[decay_model] = (table[:, 0], columns)
    return _BRANCHING_CACHE[decay_model]


def _width_table(name, decay_model=DEFAULT_DECAY_MODEL):
    """(m_a, Gamma/(1/f)^2) for one canonical width column."""
    mass_grid, columns = _load_width_tables(decay_model)
    return mass_grid, columns[name]


def _table_width(name, m_a, decay_model=DEFAULT_DECAY_MODEL):
    """Interpolated Gamma/(1/f)^2 [GeV^3]; zero outside the tabulated range."""
    m_grid, g_grid = _width_table(name, decay_model)
    m_a = float(m_a)
    if m_a < m_grid[0] or m_a > m_grid[-1]:
        return 0.0
    return float(np.interp(m_a, m_grid, g_grid))


def table_mass_max(name="total", decay_model=DEFAULT_DECAY_MODEL):
    return float(_width_table(name, decay_model)[0][-1])


def table_mass_min(name="total", decay_model=DEFAULT_DECAY_MODEL):
    return float(_width_table(name, decay_model)[0][0])


def _exclusive_branching(
    channel_id, m_a, decay_model=DEFAULT_DECAY_MODEL
):
    mass_grid, columns = _load_branching_tables(decay_model)
    m_a = float(m_a)
    if m_a < mass_grid[0] or m_a > mass_grid[-1]:
        return 0.0
    return float(np.interp(m_a, mass_grid, columns[channel_id]))


# --------------------------------------------------------------------------
# a -> f fbar partial widths, total width, lifetime
# --------------------------------------------------------------------------

def _two_body_fermion_width(m_a, m_f, n_c, inv_f, c_f=1.0):
    """Gamma(a -> f fbar) [GeV] for the pseudoscalar coupling g = c_f m_f / f."""
    m_a = float(m_a)
    if m_a <= 2.0 * m_f:
        return 0.0
    beta = np.sqrt(1.0 - (2.0 * m_f / m_a) ** 2)
    return n_c * (c_f * m_f) ** 2 * m_a * inv_f ** 2 * beta / (8.0 * np.pi)


def alp_partial_widths(
    m_a, inv_f, c_f=1.0, decay_model=DEFAULT_DECAY_MODEL
):
    """Per-channel partial widths [GeV] as a dict.

    Every channel comes from the selected exact SensCalc table in the BNT
    convention. Universal ``c_f`` and ``1/f`` rescale all coefficients together.
    """
    scale = inv_f ** 2 * c_f ** 2
    return {
        "ee": _table_width("ee", m_a, decay_model) * scale,
        "mumu": _table_width("mumu", m_a, decay_model) * scale,
        "tautau": _table_width("tautau", m_a, decay_model) * scale,
        "gammagamma": _table_width("gammagamma", m_a, decay_model) * scale,
        "hadronic": _table_width("hadronic_total", m_a, decay_model) * scale,
    }


def hadronic_invisible_width(
    m_a, inv_f, c_f=1.0, decay_model=DEFAULT_DECAY_MODEL
):
    """Hadronic width outside SensCalc's charged/no-ECAL channel selection."""
    total = alp_total_width(m_a, inv_f, c_f, decay_model)
    visible_br = sum(
        _exclusive_branching(channel_id, m_a, decay_model)
        for channel_id in SENSCALC_VISIBLE_HADRONIC_CHANNEL_IDS
    )
    hadronic = alp_partial_widths(m_a, inv_f, c_f, decay_model)["hadronic"]
    return max(hadronic - visible_br * total, 0.0)


def alp_total_width(m_a, inv_f, c_f=1.0, decay_model=DEFAULT_DECAY_MODEL):
    return float(
        _table_width("total", m_a, decay_model) * inv_f ** 2 * c_f ** 2
    )


def alp_ctau(m_a, inv_f, c_f=1.0, decay_model=DEFAULT_DECAY_MODEL):
    """Lab-frame c*tau [m] of the ALP at rest-frame width Gamma_tot."""
    gamma = alp_total_width(m_a, inv_f, c_f, decay_model)
    if gamma <= 0.0:
        return np.inf
    return HBAR_C_GEV_M / gamma


def alp_branchings(
    m_a, inv_f=INV_F_REF, c_f=1.0, decay_model=DEFAULT_DECAY_MODEL
):
    """Branching ratios per channel (coupling-independent ratios)."""
    w = alp_partial_widths(m_a, inv_f, c_f, decay_model)
    tot = alp_total_width(m_a, inv_f, c_f, decay_model)
    if tot <= 0.0:
        return {k: 0.0 for k in w}
    return {k: v / tot for k, v in w.items()}


# --------------------------------------------------------------------------
# Visible decay channels for the charged-track reconstruction
# --------------------------------------------------------------------------
# These grouped channels support only the explicitly retained legacy proxy in
# templates.py. Publication templates instead sample the full exclusive table,
# decay unstable daughters with Pythia, and apply the exact three-body weights.
# The tau decay length (c*tau ~ 87 um * beta*gamma) is far below the 3 mm hit
# resolution, so each tau is treated as one charged track along its direction.
#
# (pdg-pair used for the two charged daughters, daughter mass).
VISIBLE_CHANNELS = {
    "ee": (11, M_E),
    "mumu": (13, M_MU),
    "tautau": (15, M_TAU),
    "hadronic": (211, M_PION),   # proxy: two leading charged hadrons ~ pions
}


def visible_channel_weights(
    m_a, inv_f=INV_F_REF, c_f=1.0, decay_model=DEFAULT_DECAY_MODEL
):
    """{channel: BR} over the track-producing channels (ratios coupling
    independent).  ``hadronic`` here is the *charged-visible* hadronic BR:
    total hadronic minus the all-neutral modes."""
    if alp_total_width(m_a, inv_f, c_f, decay_model) <= 0.0:
        return {}
    out = {
        channel: value
        for channel, channel_id in (
            ("ee", "channel_001"),
            ("mumu", "channel_002"),
            ("tautau", "channel_003"),
        )
        if (value := _exclusive_branching(channel_id, m_a, decay_model)) > 0.0
    }
    hadronic = sum(
        _exclusive_branching(channel_id, m_a, decay_model)
        for channel_id in SENSCALC_VISIBLE_HADRONIC_CHANNEL_IDS
    )
    if hadronic > 0.0:
        out["hadronic"] = hadronic
    return out


def visible_fraction(
    m_a, inv_f=INV_F_REF, c_f=1.0, decay_model=DEFAULT_DECAY_MODEL
):
    """Total BR into channels with >= 2 prompt charged tracks (multiplies the
    production yield in the sensitivity scan; coupling independent)."""
    return float(
        sum(visible_channel_weights(m_a, inv_f, c_f, decay_model).values())
    )


# --------------------------------------------------------------------------
# Production: BR(B -> K^(i) a) over the kaon tower (b -> s a penguin)
# --------------------------------------------------------------------------
# Amplitudes as in ALPINIST alp_production_BRs/ds_production_BRs (BSD-3;
# form factors from Boiarska et al. arXiv:1904.10447):
#
#   Gamma(B -> K_i a) = |g_bs|^2 m_bs2(K_i) |M_i(m_a)|^2
#                        * sqrt(lambda(m_B^2, m_K^2, m_a^2)) / (16 pi m_B^3)
#
# with g_bs = CBS_EFF * (1/f), m_bs2 the (m_b -+ m_s)^2 quark-density factor
# (pseudoscalar/axial daughters: -, scalar/vector/tensor: +).

M_B_QUARK = QUARK_MASSES["b"]
M_S_QUARK = QUARK_MASSES["s"]

# Kaon-tower masses (GeV): charged partner for B+, neutral partner for B0
# where they differ (ALPINIST exotic_constants m_K dict / PDG).
KAON_TOWER = {
    #  name          m(charged)  m(neutral)
    "K":            (0.493677, 0.497611),
    "K0star_700":   (0.845, 0.845),
    "K0star_1430":  (1.425, 1.425),
    "Kstar_892":    (0.89176, 0.89555),
    "Kstar_1410":   (1.414, 1.414),
    "Kstar_1680":   (1.718, 1.718),
    "K1_1270":      (1.253, 1.253),
    "K1_1400":      (1.403, 1.403),
    "K2star_1430":  (1.4273, 1.4324),
}


def _kallen(a, b, c):
    return a * a + b * b + c * c - 2.0 * (a * b + b * c + c * a)


def _lambda_kallen_m(m1, m2, m3):
    """Kallen lambda of the squared masses, lambda(m1^2, m2^2, m3^2)."""
    return _kallen(m1 * m1, m2 * m2, m3 * m3)


def g_bs(inv_f):
    """Flavour-violating b-s-a coupling [GeV^-1] at inverse decay constant
    1/f (linear rescale off the Lambda = 1 TeV evaluation of CBS_EFF)."""
    return CBS_EFF * inv_f


def _M_BP(m_B, m_K, q, name):
    """B -> K (pseudoscalar) via the scalar density, f_0 form factor."""
    f_0, m_fit = 0.33, 6.12
    return 0.5 * (m_B ** 2 - m_K ** 2) / (M_B_QUARK - M_S_QUARK) \
        * f_0 / (1.0 - q ** 2 / m_fit ** 2)


def _M_BS(m_B, m_K, q, name):
    """B -> K0* (scalar) via the pseudoscalar density."""
    if name == "K0star_700":
        f_p, a, b = 0.46, 1.6, 1.35
    else:  # K0star_1430
        f_p, a, b = 0.17, 4.4, 6.4
    return 0.5 * (m_B ** 2 - m_K ** 2 - q ** 2) / (M_B_QUARK + M_S_QUARK) \
        * f_p / (1.0 - a * q ** 2 / m_B ** 2 + b * (q ** 2 / m_B ** 2) ** 2)


def _M_BV(m_B, m_K, q, name):
    """B -> K* (vector), A_0 form factor."""
    if name == "Kstar_892":
        a_0 = 1.364 / (1.0 - q ** 2 / m_B ** 2) - 0.990 / (1.0 - q ** 2 / 36.78)
    elif name == "Kstar_1410":
        a_0 = ((1.0 - 2.0 * m_K ** 2 / (m_B ** 2 + m_K ** 2 - q ** 2)) * 0.22
               + m_K / m_B * 0.28) / (1.0 - q ** 2 / m_B ** 2)
    else:  # Kstar_1680
        a_0 = ((1.0 - 2.0 * m_K ** 2 / (m_B ** 2 + m_K ** 2 - q ** 2)) * 0.18
               + m_K / m_B * 0.24) / (1.0 - q ** 2 / m_B ** 2)
    return -0.5 * np.sqrt(max(_lambda_kallen_m(m_B, m_K, q), 0.0)) \
        / (M_B_QUARK + M_S_QUARK) * a_0


def _M_BA(m_B, m_K, q, name):
    """B -> K1 (axial vector)."""
    th = -0.593

    def v_ab(f0, a, b):
        return f0 / (1.0 - a * (q / m_B) ** 2 + b * (q / m_B) ** 4)

    if name == "K1_1270":
        v_0 = np.sin(th) * 1.31 * v_ab(0.22, 2.4, 1.78) \
            + np.cos(th) * 1.34 * v_ab(-0.45, 1.34, 0.64)
    else:  # K1_1400
        v_0 = np.cos(th) * 1.31 * v_ab(0.22, 2.4, 1.78) \
            - np.sin(th) * 1.34 * v_ab(-0.45, 1.34, 0.64)
    return 0.5 * np.sqrt(max(_lambda_kallen_m(m_B, m_K, q), 0.0)) \
        / (M_B_QUARK - M_S_QUARK) * v_0 / m_K


def _M_BT(m_B, m_K, q, name):
    """B -> K2* (tensor)."""
    a_t = 0.23 / ((1.0 - (q / m_B) ** 2)
                  * (1.0 - 1.23 * (q / m_B) ** 2 + 0.76 * (q / m_B) ** 4))
    return -0.5 * np.sqrt(1.0 / 6.0) * (_lambda_kallen_m(m_B, m_K, q) / (m_B * m_K)) \
        / (M_B_QUARK + M_S_QUARK) * a_t


# name -> (matrix element, quark-density sign: -1 => (m_b - m_s), +1 => (m_b + m_s))
_TOWER_AMPLITUDES = {
    "K":            (_M_BP, -1),
    "K0star_700":   (_M_BS, +1),
    "K0star_1430":  (_M_BS, +1),
    "Kstar_892":    (_M_BV, +1),
    "Kstar_1410":   (_M_BV, +1),
    "Kstar_1680":   (_M_BV, +1),
    "K1_1270":      (_M_BA, -1),
    "K1_1400":      (_M_BA, -1),
    "K2star_1430":  (_M_BT, +1),
}


def _parent_props(parent):
    if parent == "B+":
        return M_BPLUS, TAU_BPLUS_S, 0
    if parent == "B0":
        return M_B0, TAU_B0_S, 1
    raise ValueError(f"unknown parent {parent!r}")


def kaon_mass(kaon, parent="B+"):
    _, _, col = _parent_props(parent)
    return KAON_TOWER[kaon][col]


def width_B_to_Ka(m_a, inv_f, parent="B+", kaon="K"):
    """Gamma(B -> K_i a) [GeV] for one kaon-tower channel."""
    m_a = float(m_a)
    m_B, _, _ = _parent_props(parent)
    m_K = kaon_mass(kaon, parent)
    if m_a >= m_B - m_K:
        return 0.0
    lam = _lambda_kallen_m(m_B, m_K, m_a)
    if lam <= 0.0:
        return 0.0
    amp_fn, sign = _TOWER_AMPLITUDES[kaon]
    m_q = M_B_QUARK + sign * M_S_QUARK
    amp2 = (g_bs(inv_f) * m_q * amp_fn(m_B, m_K, m_a, kaon)) ** 2
    return amp2 * np.sqrt(lam) / (16.0 * np.pi * m_B ** 3)


def br_B_to_Ka(m_a, inv_f, parent="B+", kaon="K"):
    """Branching ratio of one B -> K_i a channel."""
    _, tau_s, _ = _parent_props(parent)
    return width_B_to_Ka(m_a, inv_f, parent, kaon) / (HBAR_GEV_S / tau_s)


def br_B_to_K_a(m_a, inv_f, parent="B+"):
    """BR of the ground-state B -> K a channel (kept for tests/back-compat)."""
    return br_B_to_Ka(m_a, inv_f, parent, "K")


def br_B_tower_total(m_a, inv_f, parent="B+"):
    """Sum of BR(B -> K_i a) over the whole kaon tower."""
    return sum(br_B_to_Ka(m_a, inv_f, parent, k) for k in KAON_TOWER)


# Production channels: (parent label, FONLL species pdg) x kaon tower.
# Fragmentation fractions reuse the HNL convention (FRAG_B).
PRODUCTION_PARENTS = (
    ("B+", 521),
    ("B0", 511),
)
