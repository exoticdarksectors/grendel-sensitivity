#!/usr/bin/env python3
"""Parametric charged-kaon -> HNL production CSV generation."""

import random
import sys
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config_mass_grid import MASS_GRID
from production.constants import (
    M_KAON, LEPTON_MASSES, FLAVOR_TO_LEPTON_PDG, SIGMA_KAON_PB,
    KAON_TSALLIS_T, KAON_TSALLIS_N, KAON_PT_MAX, KAON_RAPIDITY_SIGMA,
    KAON_E_MAX,
)
from production.decay_engine.generate_meson_csvs import (
    compute_production_br_components, _sample_hnl_from_mesons,
)
from production.hnlcalc import init_hnlcalc
from production.io import llp_csv_path, write_empty_csv, write_llp_csv
from production.paths import LLP_VECTORS_DIR

OUTPUT_BASE = LLP_VECTORS_DIR

N_POOL = 100_000
KAON_PDG = 321


def sample_kaon_4vectors(n_pool, rng):
    m_k = M_KAON

    n_grid = 4000
    pt_edges = np.linspace(0.0, KAON_PT_MAX, n_grid + 1)
    pt_centers = 0.5 * (pt_edges[:-1] + pt_edges[1:])
    mt = np.sqrt(pt_centers**2 + m_k**2)
    pdf = pt_centers * (1.0 + (mt - m_k) / (KAON_TSALLIS_N * KAON_TSALLIS_T)) ** (-KAON_TSALLIS_N)
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]

    u = rng.random(n_pool)
    idx = np.searchsorted(cdf, u)
    idx = np.clip(idx, 0, n_grid - 1)
    pt = rng.uniform(pt_edges[idx], pt_edges[idx + 1])

    y = rng.normal(0.0, KAON_RAPIDITY_SIGMA, n_pool)
    mt = np.sqrt(pt * pt + m_k * m_k)
    bad = mt * np.cosh(y) >= KAON_E_MAX
    while bad.any():
        y[bad] = rng.normal(0.0, KAON_RAPIDITY_SIGMA, int(bad.sum()))
        bad = mt * np.cosh(y) >= KAON_E_MAX
    phi = rng.uniform(0.0, 2.0 * np.pi, n_pool)

    mt_samp = np.sqrt(pt**2 + m_k**2)
    pz = mt_samp * np.sinh(y)
    E = mt_samp * np.cosh(y)
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)

    return {
        'E': E, 'px': px, 'py': py, 'pz': pz,
        'species_pdg': np.full(n_pool, KAON_PDG, dtype=int),
        'pt': pt, 'y': y, 'phi': phi,
    }


def process_flavor(flavor, pool, masses, rng):
    hnl = init_hnlcalc(flavor)
    lepton_pdg = FLAVOR_TO_LEPTON_PDG[flavor]
    m_lepton = LEPTON_MASSES[flavor]

    n_pool = len(pool['E'])
    m_parent = M_KAON

    for m_N in masses:
        csv_path = llp_csv_path(flavor, "Kmeson", m_N, base=OUTPUT_BASE)

        if m_N >= m_parent - m_lepton:
            write_empty_csv(csv_path)
            continue

        _, br_3body_channels, br = compute_production_br_components(
            hnl, KAON_PDG, lepton_pdg, m_N
        )
        if br <= 0:
            write_empty_csv(csv_path)
            continue

        hnl_4v = _sample_hnl_from_mesons(
            pool['E'], pool['px'], pool['py'], pool['pz'],
            m_parent, m_lepton, m_N,
            br_3body_channels, br, hnl, rng,
        )
        w = SIGMA_KAON_PB * br / n_pool
        weights = np.full(n_pool, w)
        write_llp_csv(csv_path, weights, hnl_4v[:, 0], hnl_4v[:, 1], hnl_4v[:, 2], hnl_4v[:, 3])
        print(f"    {csv_path.name}: {n_pool} events, BR={br:.3e}, w_sum={weights.sum():.3e}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Generate kaon -> HNL CSVs (K+ -> l N)")
    parser.add_argument("--flavor", choices=["Ue", "Umu", "Utau"], nargs="+",
                        default=["Ue", "Umu", "Utau"])
    parser.add_argument("--n-pool", type=int, default=N_POOL)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--masses", type=float, nargs="+", default=None,
                        help="Custom mass list (default: full grid)")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    random.seed(args.seed)
    masses = args.masses if args.masses else MASS_GRID

    print(f"Sampling charged-kaon pool ({args.n_pool} events)...")
    pool = sample_kaon_4vectors(args.n_pool, rng)

    for flavor in args.flavor:
        print(f"\n{'='*60}\nFlavor: {flavor}\n{'='*60}")
        process_flavor(flavor, pool, masses, rng)

    print("\nDone.")


if __name__ == "__main__":
    main()
