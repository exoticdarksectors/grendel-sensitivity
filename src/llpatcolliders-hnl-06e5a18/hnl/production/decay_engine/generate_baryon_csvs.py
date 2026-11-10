#!/usr/bin/env python3
"""b-baryon -> HNL production CSV generation (Lambda_b -> Lambda_c l N).

The dominant baryonic HNL-production channel. The closure remainder of the
LHCb bottom-fragmentation ratios is represented as Lambda_b-like production;
Xi_b/Omega_b are not separately modeled. The Lambda_b pT-y shape reuses the
bottom FONLL grid evaluated with the Lambda_b mass (the same approximation used
for Bc). The differential rate is parametrized in (m12sq, q2) and sampled with
``decay_3body_weighted_dq2dm122``; d3 is the HNL.
"""

import random
import sys
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config_mass_grid import MASS_GRID
from production.constants import (
    M_LAMBDA_B, M_LAMBDA_C, FRAG_LAMBDA_B,
    LEPTON_MASSES, FLAVOR_TO_LEPTON_PDG,
)
from production.fonll.fonll_parser import get_sigma_total
from production.fonll.meson_sampler import (
    sample_meson_4vectors, meson_4vec_from_kinematics,
)
from production.decay_engine.kinematics import decay_3body_weighted_dq2dm122
from production.hnlcalc import init_hnlcalc
from production.io import llp_csv_path, write_empty_csv, write_llp_csv
from production.paths import LLP_VECTORS_DIR

OUTPUT_BASE = LLP_VECTORS_DIR

N_POOL = 100_000

LAMBDA_B_PDG = 5122
LAMBDA_C_PDG = 4122


def build_lambda_b_pool(n_pool, rng):
    """Sample a Lambda_b lab-frame 4-vector pool from the bottom FONLL shape.

    Returns ``(pool_dict, sigma_bottom)`` where ``pool_dict`` has the same
    layout as the meson pools (``E, px, py, pz, pt, y, phi``) but the four-vector
    is built with the Lambda_b mass. ``sigma_bottom`` is the inclusive bottom
    FONLL cross section in pb (b only; the b+bbar factor of 2 is applied in the
    per-mass weight, matching the meson channels).
    """
    # Reuse the bottom FONLL (pT, y) sampling; rebuild the 4-vector with M_LAMBDA_B.
    raw = sample_meson_4vectors(n_pool, "bottom", rng=rng)
    v = meson_4vec_from_kinematics(raw['pt'], raw['y'], raw['phi'], M_LAMBDA_B)
    pool = {
        'E': v['E'], 'px': v['px'], 'py': v['py'], 'pz': v['pz'],
        'species_pdg': np.full(n_pool, LAMBDA_B_PDG, dtype=int),
        'pt': raw['pt'], 'y': raw['y'], 'phi': raw['phi'],
    }
    sigma_bottom = get_sigma_total("bottom")
    return pool, sigma_bottom


def process_flavor(flavor, pool, sigma_bottom, masses, rng):
    hnl = init_hnlcalc(flavor)
    lepton_pdg = FLAVOR_TO_LEPTON_PDG[flavor]
    m_lepton = LEPTON_MASSES[flavor]

    n_pool = len(pool['E'])
    m_parent = M_LAMBDA_B
    m_daughter = M_LAMBDA_C
    threshold = M_LAMBDA_B - M_LAMBDA_C - m_lepton

    for m_N in masses:
        csv_path = llp_csv_path(flavor, "Bbaryon", m_N, base=OUTPUT_BASE)

        if m_N >= threshold:
            write_empty_csv(csv_path)
            continue

        dbr = hnl.get_3body_dbr_baryon(
            str(LAMBDA_B_PDG), str(LAMBDA_C_PDG), str(lepton_pdg)
        )
        br = hnl.integrate_3body_br(
            dbr, m_N, M_LAMBDA_B, M_LAMBDA_C, m_lepton,
            coupling=1.0, nsample=500, integration='dq2dm122',
        )
        if br is None or np.isnan(br) or br <= 0:
            write_empty_csv(csv_path)
            continue

        # d1 = Lambda_c, d2 = lepton, d3 = N (the HNL).
        _, _, hnl_4v = decay_3body_weighted_dq2dm122(
            pool['E'], pool['px'], pool['py'], pool['pz'],
            m_parent, m_daughter, m_lepton, m_N,
            dbr_expr=dbr, coupling=1.0, rng=rng,
        )

        # b + bbar convention (factor 2): Lambda_b from b, Lambda_b-bar from bbar.
        w = 2.0 * sigma_bottom * FRAG_LAMBDA_B * br / n_pool
        weights = np.full(n_pool, w)
        write_llp_csv(
            csv_path, weights,
            hnl_4v[:, 0], hnl_4v[:, 1], hnl_4v[:, 2], hnl_4v[:, 3],
        )
        print(f"    {csv_path.name}: {n_pool} events, BR={br:.3e}, "
              f"w_sum={weights.sum():.3e}")


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Generate b-baryon -> HNL CSVs (Lambda_b -> Lambda_c l N)"
    )
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

    print(f"Sampling Lambda_b pool ({args.n_pool} events)...")
    pool, sigma_bottom = build_lambda_b_pool(args.n_pool, rng)
    print(f"  sigma_bottom(FONLL) = {sigma_bottom:.3e} pb")

    for flavor in args.flavor:
        print(f"\n{'='*60}\nFlavor: {flavor}\n{'='*60}")
        process_flavor(flavor, pool, sigma_bottom, masses, rng)

    print("\nDone.")


if __name__ == "__main__":
    main()
