"""BC5 production: the S four-vectors of every channel, one CSV per mass."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ...io.vectors import csv_column_count, format_mass_for_filename, write_llp_csv
from ...production.decay_engine.kinematics import decay_2body
from ...production.fonll.fonll_parser import get_sigma_total
from ...production.fonll.meson_sampler import meson_4vec_from_kinematics, sample_meson_4vectors
from . import model, production
from .higgs_pool import sample_higgs_4vectors

MASS_GRID_BC5 = sorted({round(x, 3) for x in production.MASS_GRID} | {
    4.8, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0, 15.0, 16.0, 18.0, 20.0,
    22.5, 25.0, 27.5, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 62.0})
N_HIGGS_DEFAULT = 200_000
QUARTIC_POOL_FRACTION = 0.25
# The production modes, in row order: b -> X_s S through the mixing angle (the BC4 channel, weight
# scaling with sin^2 theta), on-shell h -> S S, and b -> X_s S S plus B_s -> S S through the off-shell
# Higgs (the last two fixed by BR(h -> SS)). Each row's mode is its ``channel`` column.
CHANNELS = ("mixing", "hSS", "BSS")
CHANNEL_ID = {name: index for index, name in enumerate(CHANNELS)}
N_COLUMNS = 7


def _empty():
    return (np.empty(0),) * 5


def _both(s1, s2, w):
    """Stack the two scalars of a pair, each carrying the pair's weight."""
    p4 = np.concatenate([s1, s2], axis=0)
    return np.concatenate([w, w]), p4[:, 0], p4[:, 1], p4[:, 2], p4[:, 3]


def generate_hSS_4vectors(m_S, higgs, rng, br_hss=model.BR_HSS_BC5):
    """S four-vectors of h -> S S from a Higgs sample ``higgs`` (the dict of
    ``higgs_pool.sample_higgs_4vectors``): weight sigma_h * BR / n per scalar, two per Higgs."""
    n = len(higgs["E"])
    if m_S >= model.M_S_MAX_HSS or br_hss <= 0.0 or n == 0:
        return _empty()
    s1, s2 = decay_2body(higgs["E"], higgs["px"], higgs["py"], higgs["pz"],
                         model.M_HIGGS, m_S, m_S, rng=rng)
    return _both(s1, s2, np.full(n, higgs["sigma_pb"] * br_hss / n))


def _pool_subset(pool, fraction):
    """The first ``fraction`` of the (i.i.d.) pool entries."""
    n = max(1, int(round(len(pool["pt"]) * fraction)))
    return {key: (value[:n] if isinstance(value, np.ndarray) and value.shape[:1] == pool["pt"].shape[:1]
                  else value) for key, value in pool.items()}


def generate_quartic_b_4vectors(m_S, pool, rng, alpha, sigma_bottom, pool_fraction=QUARTIC_POOL_FRACTION):
    """S four-vectors of b-hadron -> X_s S S and B_s -> S S at the quartic coupling ``alpha`` (GeV),
    from the first ``pool_fraction`` of the shared b-hadron ``pool``."""
    if alpha <= 0.0:
        return _empty()
    pool = _pool_subset(pool, pool_fraction)
    n_each = len(pool["pt"])
    sampling_weight = pool.get("sampling_weight", np.ones(n_each))
    parts = []
    for parent, m_B, m_recoil, frag in production.B_SPECIES.values():
        v = meson_4vec_from_kinematics(pool["pt"], pool["y"], pool["phi"], m_B)
        if m_S < 0.5 * (m_B - m_recoil):
            br = float(model.br_b_to_Xs_SS(m_S, parent=parent, alpha=alpha))
            q2 = model.sample_q2_B_to_K_SS(n_each, m_S, rng, m_B, m_recoil)
            if br > 0.0 and len(q2) == n_each:
                q = np.sqrt(q2)
                _, ss = decay_2body(v["E"], v["px"], v["py"], v["pz"], m_B, m_recoil, q, rng=rng)
                s1, s2 = decay_2body(ss[:, 0], ss[:, 1], ss[:, 2], ss[:, 3], q, m_S, m_S, rng=rng)
                w = 2.0 * sigma_bottom * frag * br / n_each * sampling_weight
                parts.append(_both(s1, s2, w))
        if parent == "Bs" and m_S < 0.5 * m_B:
            br = float(model.br_Bs_to_SS(m_S, alpha=alpha))
            if br > 0.0:
                s1, s2 = decay_2body(v["E"], v["px"], v["py"], v["pz"], m_B, m_S, m_S, rng=rng)
                w = 2.0 * sigma_bottom * frag * br / n_each * sampling_weight
                parts.append(_both(s1, s2, w))
    if not parts:
        return _empty()
    return tuple(np.concatenate([part[i] for part in parts]) for i in range(5))


def generate_bc5_4vectors(m_S, rngs, *, pool, higgs, sigma_bottom, br_hss=model.BR_HSS_BC5,
                          quartic_b=True):
    """Every BC5 row for ``m_S``: ``(weights, E, px, py, pz, coupling_power, channel)`` plus the row
    count per channel."""
    rng_mix, rng_q, rng_h = rngs
    n_pool = len(pool["pt"])
    parts = {"mixing": production.generate_scalar_4vectors(
        m_S, n_pool, rng_mix, sigma_bottom=sigma_bottom, pool=pool)}
    alpha = model.alpha_quartic(m_S, br_hss)
    parts["hSS"] = (generate_hSS_4vectors(m_S, higgs, rng_h, br_hss)
                    if higgs is not None else _empty())
    parts["BSS"] = (generate_quartic_b_4vectors(m_S, pool, rng_q, alpha, sigma_bottom)
                    if quartic_b else _empty())
    counts = {name: len(parts[name][0]) for name in CHANNELS}
    power = np.concatenate([np.ones(counts["mixing"]), np.zeros(counts["hSS"]), np.zeros(counts["BSS"])])
    channel = np.concatenate([np.full(counts[name], CHANNEL_ID[name], float) for name in CHANNELS])
    arrays = tuple(np.concatenate([parts[name][i] for name in CHANNELS]) for i in range(5))
    return (*arrays, power, channel), counts


def write_bc5_csv(m_S, out_dir, rngs, **kwargs):
    """Generate and write ``<out_dir>/mS_<label>.csv`` (``N_COLUMNS`` columns)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"mS_{format_mass_for_filename(m_S)}.csv"
    (w, E, px, py, pz, power, channel), counts = generate_bc5_4vectors(m_S, rngs, **kwargs)
    if len(w) == 0:
        path.write_text("")
    else:
        write_llp_csv(path, w, E, px, py, pz, coupling_power=power, channel=channel)
    return path, counts


def is_current_format(path) -> bool:
    """Whether an existing BC5 four-vector CSV has every column of the current format (the channel
    column arrived after v2); an empty CSV, a mass with no rows, has nothing to add."""
    n = csv_column_count(path)
    return n == 0 or n >= N_COLUMNS


def make_streams(seed):
    """The three generators of a production run (mixing, quartic B, Higgs)."""
    return (np.random.default_rng(seed), np.random.default_rng(seed + 1),
            np.random.default_rng(seed + 2))


def shared_pools(rngs, n_pool, n_higgs, *, br_hss=model.BR_HSS_BC5,
                 high_pt_tilt_scale=None, nominal_mixture_fraction=0.5):
    """The b-hadron pool (mixing stream, as BC4) and the Higgs sample (Higgs stream) every mass point
    decays; ``(pool, higgs, sigma_bottom)``."""
    rng_mix, _, rng_h = rngs
    pool = sample_meson_4vectors(n_pool, "bottom", rng=rng_mix,
                                 high_pt_tilt_scale=high_pt_tilt_scale,
                                 nominal_mixture_fraction=nominal_mixture_fraction)
    higgs = sample_higgs_4vectors(n_higgs, rng_h) if (br_hss > 0.0 and n_higgs > 0) else None
    return pool, higgs, get_sigma_total("bottom")
