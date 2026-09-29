r"""BC10 fermiophilic-ALP model layer (couplings -> production, lifetime, BRs)."""

from __future__ import annotations

import json

import numpy as np
from particle import Particle

from ...io.tables import open_table, read_table_text
from .decay_models import DEFAULT_DECAY_MODEL, decay_data_dir

HBAR_C_GEV_M = 1.973269804e-16
HBAR_GEV_S = 6.582119569e-25

M_E = Particle.from_pdgid(11).mass * 1e-3
M_MU = Particle.from_pdgid(13).mass * 1e-3
M_TAU = Particle.from_pdgid(15).mass * 1e-3
M_PION = Particle.from_pdgid(211).mass * 1e-3

QUARK_MASSES = {
    "u": 0.00216, "d": 0.00467, "s": 0.0934, "c": 1.27, "b": 4.18,
}

M_BPLUS = Particle.from_pdgid(521).mass * 1e-3
M_B0 = Particle.from_pdgid(511).mass * 1e-3
TAU_BPLUS_S = 1.638e-12
TAU_B0_S = 1.519e-12

CBS_EFF = 3.518383e-4
INV_F_REF = 1.0e-3

LIGHT_MESON_RESONANCE_WINDOWS = (
    (0.125, 0.140, "pi0"),
    (0.538, 0.555, "eta"),
    (0.940, 0.974, "eta-prime"),
)

SENSCALC_2501_DATA_DIR = decay_data_dir(DEFAULT_DECAY_MODEL)

_WIDTH_CACHE: dict[str, tuple[np.ndarray, dict[str, np.ndarray]]] = {}
_BRANCHING_CACHE: dict[str, tuple[np.ndarray, dict[str, np.ndarray]]] = {}

SENSCALC_INVISIBLE_CHANNEL_IDS = {
    "channel_004",
    "channel_014",
    "channel_021",
}
SENSCALC_DUPLICATE_CHANNEL_IDS = {
    "channel_024",
    "channel_029",
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
    """Per-channel partial widths [GeV] as a dict."""
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
    """Proper decay length c*tau = hbar c / Gamma_tot [m] of the ALP."""
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


VISIBLE_CHANNELS = {
    "ee": (11, M_E),
    "mumu": (13, M_MU),
    "tautau": (15, M_TAU),
    "hadronic": (211, M_PION),
}


def visible_channel_weights(
    m_a, inv_f=INV_F_REF, c_f=1.0, decay_model=DEFAULT_DECAY_MODEL
):
    """{channel: BR} over the track-producing channels (ratios coupling independent)."""
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
    """Total BR of SensCalc's no-ECAL visible selection (multiplies the production yield only for
    template bundles without the full branching, i.e."""
    return float(
        sum(visible_channel_weights(m_a, inv_f, c_f, decay_model).values())
    )


M_B_QUARK = QUARK_MASSES["b"]
M_S_QUARK = QUARK_MASSES["s"]

KAON_TOWER = {
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
    """Flavour-violating b-s-a coupling [GeV^-1] at inverse decay constant 1/f (linear rescale off the
    Lambda = 1 TeV evaluation of CBS_EFF)."""
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
    else:
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
    else:
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
    else:
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


PRODUCTION_PARENTS = (
    ("B+", 521),
    ("B0", 511),
)
