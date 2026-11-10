#!/usr/bin/env python3
"""B/D/Bc -> HNL production CSV generation."""

import random
import sys
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config_mass_grid import MASS_GRID
from production.constants import (
    MESON_MASSES, LEPTON_MASSES, FLAVOR_TO_LEPTON_PDG,
    QUARK_MESON_MAP, SIGMA_BC_PB,
)
from production.fonll.fonll_parser import get_sigma_total
from production.fonll.meson_sampler import (
    sample_meson_4vectors, meson_4vec_from_kinematics,
)
from production.decay_engine.kinematics import (
    decay_2body, decay_3body_weighted_dq2dE,
)
from production.hnlcalc import init_hnlcalc
from production.io import llp_csv_path, write_empty_csv, write_llp_csv
from production.paths import LLP_VECTORS_DIR

OUTPUT_BASE = LLP_VECTORS_DIR

N_POOL = 100_000

CHANNEL_LABELS = {
    "bottom": "Bmeson",
    "charm": "Dmeson",
    "bc": "Bc",
}

CHARGED_PSEUDOSCALARS = {211, 321, 411, 431, 521, 541}

THREEBODY_CHANNELS = {
    321: [
        (111, "pseudo"),   # K+ -> pi0 l+ N
    ],
    421: [  # D0
        (-321, "pseudo"),   # D0 → K- ℓ+ N
        (-323, "vector"),   # D0 → K*- ℓ+ N
        (-211, "pseudo"),   # D0 → π- ℓ+ N
        (-213, "vector"),   # D0 → ρ- ℓ+ N
    ],
    411: [  # D+
        (-311, "pseudo"),   # D+ → K̄0 ℓ+ N
        (-313, "vector"),   # D+ → K̄*0 ℓ+ N
        (111, "pseudo"),    # D+ → π0 ℓ+ N
        (113, "vector"),    # D+ → ρ0 ℓ+ N
        (221, "pseudo"),    # D+ → η ℓ+ N
        (331, "pseudo"),    # D+ → η' ℓ+ N
    ],
    431: [  # Ds
        (221, "pseudo"),    # Ds → η ℓ+ N
        (331, "pseudo"),    # Ds → η' ℓ+ N
        (311, "pseudo"),    # Ds → K0 ℓ+ N
        (313, "vector"),    # Ds → K*0 ℓ+ N
        (333, "vector"),    # Ds → φ ℓ+ N
    ],
    521: [  # B+
        (-421, "pseudo"),   # B+ → D̄0 ℓ+ N
        (-423, "vector"),   # B+ → D̄*0 ℓ+ N
        (111, "pseudo"),    # B+ → π0 ℓ+ N
        (113, "vector"),    # B+ → ρ0 ℓ+ N
        (221, "pseudo"),    # B+ → η ℓ+ N
        (223, "vector"),    # B+ → ω ℓ+ N
        (331, "pseudo"),    # B+ → η' ℓ+ N
    ],
    511: [  # B0
        (-411, "pseudo"),   # B0 → D- ℓ+ N
        (-413, "vector"),   # B0 → D*- ℓ+ N
        (-211, "pseudo"),   # B0 → π- ℓ+ N
        (-213, "vector"),   # B0 → ρ- ℓ+ N
    ],
    531: [  # Bs
        (-431, "pseudo"),   # Bs → Ds- ℓ+ N
        (-433, "vector"),   # Bs → Ds*- ℓ+ N
        (-321, "pseudo"),   # Bs → K- ℓ+ N
        (-323, "vector"),   # Bs → K*- ℓ+ N
    ],
    541: [  # Bc
        (511, "pseudo"),    # Bc → B0 ℓ+ N
        (531, "pseudo"),    # Bc → Bs ℓ+ N
        (513, "vector"),    # Bc → B*0 ℓ+ N
        (533, "vector"),    # Bc → Bs* ℓ+ N
        (421, "pseudo"),    # Bc → D0 ℓ+ N
        (441, "pseudo"),    # Bc → ηc ℓ+ N
        (423, "vector"),    # Bc → D*0 ℓ+ N
        (443, "vector"),    # Bc → J/ψ ℓ+ N
    ],
}


def _eval_2body_br(hnl, parent_pdg, lepton_pdg, m_N):
    if abs(parent_pdg) not in CHARGED_PSEUDOSCALARS:
        return 0.0
    sign = "-" if parent_pdg > 0 else ""
    pid_lep = f"{sign}{abs(lepton_pdg)}"
    br_expr = hnl.get_2body_br(str(parent_pdg), pid_lep)
    br_val = eval(br_expr, {"np": np, "__builtins__": {}},
                  {"mass": m_N, "coupling": 1.0})
    if np.isnan(br_val) or br_val < 0:
        return 0.0
    return float(br_val)


def _eval_3body_br(hnl, parent_pdg, daughter_pdg, lepton_pdg, m_N, ch_type):
    m_parent = MESON_MASSES.get(abs(parent_pdg), hnl.masses(parent_pdg))
    m_daughter = hnl.masses(daughter_pdg)
    m_lepton = hnl.masses(lepton_pdg)

    if m_N >= m_parent - m_daughter - m_lepton:
        return 0.0, None

    sign_lep = "-" if parent_pdg > 0 else ""
    pid_lep_str = f"{sign_lep}{abs(lepton_pdg)}"

    if ch_type == "pseudo":
        dbr = hnl.get_3body_dbr_pseudoscalar(str(parent_pdg), str(daughter_pdg), pid_lep_str)
    elif ch_type == "vector":
        dbr = hnl.get_3body_dbr_vector(str(parent_pdg), str(daughter_pdg), pid_lep_str)
    else:
        return 0.0, None

    br_val = hnl.integrate_3body_br(
        dbr, m_N, m_parent, m_daughter, m_lepton,
        coupling=1.0, nsample=500,
    )
    if br_val is None or np.isnan(br_val) or br_val < 0:
        return 0.0, None
    return float(br_val), dbr


def compute_total_production_br(hnl, parent_pdg, lepton_pdg, m_N):
    _, _, br_total = compute_production_br_components(hnl, parent_pdg, lepton_pdg, m_N)
    return br_total


def compute_production_br_components(hnl, parent_pdg, lepton_pdg, m_N):
    m_parent = MESON_MASSES.get(abs(parent_pdg), hnl.masses(parent_pdg))
    m_lepton = hnl.masses(lepton_pdg)

    br_2body = 0.0
    if m_N < m_parent - m_lepton:
        br_2body = _eval_2body_br(hnl, parent_pdg, lepton_pdg, m_N)

    br_3body_channels = []
    abs_pdg = abs(parent_pdg)
    if abs_pdg in THREEBODY_CHANNELS:
        for daughter_pdg, ch_type in THREEBODY_CHANNELS[abs_pdg]:
            d_pdg = -daughter_pdg if parent_pdg < 0 else daughter_pdg
            br_3b, dbr = _eval_3body_br(hnl, parent_pdg, d_pdg, lepton_pdg, m_N, ch_type)
            if br_3b > 0 and dbr is not None:
                br_3body_channels.append((d_pdg, dbr, br_3b))

    br_3body_total = sum(ch[2] for ch in br_3body_channels)
    br_total = br_2body + br_3body_total
    return br_2body, br_3body_channels, br_total


def _sample_hnl_from_mesons(parent_E, parent_px, parent_py, parent_pz, m_parent, m_lepton, m_N,
                            br_3body_channels, br_total, hnl, rng):
    n_events = len(parent_E)
    hnl_4v = np.empty((n_events, 4))

    br_3body_total = sum(ch[2] for ch in br_3body_channels)
    use_3body = np.zeros(n_events, dtype=bool)
    if br_3body_total > 0:
        p_3body = br_3body_total / br_total
        use_3body = rng.random(n_events) < p_3body

    use_2body = ~use_3body
    if use_2body.any():
        _, hnl_2b = decay_2body(
            parent_E[use_2body], parent_px[use_2body], parent_py[use_2body], parent_pz[use_2body],
            m_parent, m_lepton, m_N, rng=rng,
        )
        hnl_4v[use_2body] = hnl_2b

    if use_3body.any():
        br_arr = np.array([ch[2] for ch in br_3body_channels], dtype=float)
        prob_arr = br_arr / br_arr.sum()
        channel_idx = rng.choice(len(br_3body_channels), size=use_3body.sum(), p=prob_arr)
        evt_idx = np.where(use_3body)[0]
        for idx_ch in np.unique(channel_idx):
            daughter_pdg, dbr_expr, _ = br_3body_channels[int(idx_ch)]
            sel = evt_idx[channel_idx == idx_ch]
            m_daughter = hnl.masses(daughter_pdg)
            _, _, hnl_3b = decay_3body_weighted_dq2dE(
                parent_E[sel], parent_px[sel], parent_py[sel], parent_pz[sel],
                m_parent, m_daughter, m_lepton, m_N,
                dbr_expr=dbr_expr, coupling=1.0, rng=rng,
            )
            hnl_4v[sel] = hnl_3b

    return hnl_4v


def generate_pool(quark, n_pool, rng):
    print(f"  Generating {quark} meson pool ({n_pool} events)...")
    pool = sample_meson_4vectors(n_pool, quark, rng=rng)
    sigma = get_sigma_total(quark)
    print(f"  σ_FONLL({quark}) = {sigma:.3e} pb")
    return pool, sigma


def process_channel(flavor, quark, pool, sigma_fonll, masses, rng):
    hnl = init_hnlcalc(flavor)
    lepton_pdg = FLAVOR_TO_LEPTON_PDG[flavor]
    m_lepton = LEPTON_MASSES[flavor]
    channel_label = CHANNEL_LABELS[quark]

    n_pool = len(pool['E'])

    for m_N in masses:
        csv_path = llp_csv_path(flavor, channel_label, m_N, base=OUTPUT_BASE)

        all_weights = []
        all_E = []
        all_px = []
        all_py = []
        all_pz = []

        if quark == "bc":
            parent_pdg = 541
            m_parent = MESON_MASSES[541]
            if m_N >= m_parent - m_lepton:
                write_empty_csv(csv_path)
                continue

            _, br_3body_channels, br = compute_production_br_components(
                hnl, parent_pdg, lepton_pdg, m_N
            )
            if br <= 0:
                write_empty_csv(csv_path)
                continue

            w = SIGMA_BC_PB * br / n_pool

            mask = np.ones(n_pool, dtype=bool)
            hnl_4v = _sample_hnl_from_mesons(
                pool['E'][mask], pool['px'][mask], pool['py'][mask], pool['pz'][mask],
                m_parent, m_lepton, m_N,
                br_3body_channels, br, hnl, rng,
            )
            all_weights.append(np.full(mask.sum(), w))
            all_E.append(hnl_4v[:, 0])
            all_px.append(hnl_4v[:, 1])
            all_py.append(hnl_4v[:, 2])
            all_pz.append(hnl_4v[:, 3])

        else:
            n_pool = len(pool['pt'])
            open_species = []
            for species_pdg, frag in QUARK_MESON_MAP[quark]:
                m_parent = MESON_MASSES[species_pdg]
                if m_N >= m_parent - m_lepton or frag <= 0:
                    continue
                _, br_3body_channels, br = compute_production_br_components(
                    hnl, int(species_pdg), lepton_pdg, m_N)
                if br <= 0:
                    continue
                open_species.append((species_pdg, m_parent, frag, br, br_3body_channels))

            if open_species:
                n_each = max(1, n_pool // len(open_species))
                for species_pdg, m_parent, frag, br, br_3body_channels in open_species:
                    idx = rng.integers(0, n_pool, size=n_each)
                    v = meson_4vec_from_kinematics(
                        pool['pt'][idx], pool['y'][idx], pool['phi'][idx], m_parent)
                    w = 2.0 * sigma_fonll * frag * br / n_each
                    hnl_4v = _sample_hnl_from_mesons(
                        v['E'], v['px'], v['py'], v['pz'],
                        m_parent, m_lepton, m_N, br_3body_channels, br, hnl, rng)
                    all_weights.append(np.full(n_each, w))
                    all_E.append(hnl_4v[:, 0]); all_px.append(hnl_4v[:, 1])
                    all_py.append(hnl_4v[:, 2]); all_pz.append(hnl_4v[:, 3])

        if all_weights:
            weights = np.concatenate(all_weights)
            E = np.concatenate(all_E)
            px = np.concatenate(all_px)
            py = np.concatenate(all_py)
            pz = np.concatenate(all_pz)
            write_llp_csv(csv_path, weights, E, px, py, pz)
            print(f"    {csv_path.name}: {len(weights)} events, w_sum={weights.sum():.3e}")
        else:
            write_empty_csv(csv_path)
            print(f"    {csv_path.name}: 0 events (below threshold)")


def generate_bc_pool(n_pool, rng):
    return sample_meson_4vectors(n_pool, "bottom", rng=rng, force_species=541)


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Generate meson → HNL 4-vector CSVs")
    parser.add_argument("--flavor", choices=["Ue", "Umu", "Utau"], nargs="+",
                        default=["Ue", "Umu", "Utau"])
    parser.add_argument("--channel", choices=["Bmeson", "Dmeson", "Bc", "all"],
                        default="all")
    parser.add_argument("--n-pool", type=int, default=N_POOL)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--masses", type=float, nargs="+", default=None,
                        help="Custom mass list (default: full grid)")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    random.seed(args.seed)
    masses = args.masses if args.masses else MASS_GRID

    if args.channel == "all":
        channels = ["bottom", "charm", "bc"]
    else:
        label_to_quark = {"Bmeson": "bottom", "Dmeson": "charm", "Bc": "bc"}
        channels = [label_to_quark[args.channel]]

    pools = {}
    sigmas = {}
    for ch in channels:
        if ch == "bc":
            pools["bc"] = generate_bc_pool(args.n_pool, rng)
            sigmas["bc"] = 0.0
        else:
            pools[ch], sigmas[ch] = generate_pool(ch, args.n_pool, rng)

    for flavor in args.flavor:
        print(f"\n{'='*60}")
        print(f"Flavor: {flavor}")
        print(f"{'='*60}")
        for ch in channels:
            print(f"\n  Channel: {CHANNEL_LABELS[ch]}")
            process_channel(flavor, ch, pools[ch], sigmas.get(ch, 0.0), masses, rng)

    print("\nDone.")


if __name__ == "__main__":
    main()
