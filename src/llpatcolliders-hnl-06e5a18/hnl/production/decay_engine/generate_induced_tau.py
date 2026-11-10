#!/usr/bin/env python3
"""Heavy-hadron -> tau -> HNL production CSV generation.

The induced-tau pool is built from 11 selected heavy-hadron decays that yield a
real tau:

- two-body leptonic modes ``Ds+/D+/B+/Bc+ -> tau nu`` (fixed longitudinal
  analyzing-power approximation);
- three-body semitauonic modes ``B/Bs -> D(*)/Ds(*) tau nu`` and the baryonic
  ``Lambda_b -> Lambda_c tau nu`` (the dominant secondary tau sources).

Each source contributes a fixed block of taus whose per-event weight reproduces
its cross section, so the combined pool stays bounded at ~n_pool rows. Every
block also carries its longitudinal analyzing power, and the tau decay is run
per-polarization-block so the semitauonic/baryonic sources (whose true tau
polarization is Dalitz-dependent) are not forced to the two-body helicity value.
"""

import random
import sys
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config_mass_grid import MASS_GRID
from production.constants import (
    M_TAU, M_DS, M_DPLUS, M_BPLUS, M_B0, M_BS, M_BC,
    M_D0, M_DSTAR0, M_DSTARP, M_DSSTAR, M_LAMBDA_B, M_LAMBDA_C,
    FRAG_B, FRAG_C, FRAG_LAMBDA_B, SIGMA_BC_PB,
    BR_DS_TAUNU, BR_DPLUS_TAUNU, BR_BPLUS_TAUNU, BR_BC_TAUNU,
    BR_BP_D0_TAUNU, BR_BP_DSTAR0_TAUNU,
    BR_B0_DP_TAUNU, BR_B0_DSTARP_TAUNU,
    BR_BS_DS_TAUNU, BR_BS_DSSTAR_TAUNU, BR_LB_LC_TAUNU,
)
from production.fonll.fonll_parser import get_sigma_total
from production.fonll.meson_sampler import (
    sample_meson_4vectors, meson_4vec_from_kinematics,
)
from production.decay_engine.kinematics import (
    decay_2body, decay_3body_weighted_dq2dE, decay_3body_weighted_dq2dm122,
)
from production.decay_engine.tau_decay import (
    init_hnlcalc, compute_tau_production_br_components, sample_hnl_from_tau,
    TAU_MESON_2BODY_ASYMMETRY,
)
from production.io import llp_csv_path, write_empty_csv, write_csv_matrix
from production.paths import LLP_VECTORS_DIR

OUTPUT_BASE = LLP_VECTORS_DIR

N_POOL = 100_000

# (parent_pdg, m_parent, daughter_pdg, m_daughter, kind, quark, frag, BR, asym)
#   kind   : '2body'           parent -> tau nu (fixed longitudinal asymmetry)
#            'pseudo'/'vector'  parent -> P/V tau nu (semitauonic)
#            'baryon'           Lambda_b -> Lambda_c tau nu (sampled via dq2dm122)
#   quark  : 'charm'/'bottom' pick the FONLL shape and the 2*sigma*frag norm;
#            'bc' uses SIGMA_BC_PB directly (both charges, no factor 2, no frag).
#   asym   : longitudinal asymmetry applied to the N in the 2-body
#            tau->meson N decay. The 2-body leptonic sources produce
#            helicity-suppressed "wrong-helicity" taus, so the N is emitted
#            backward (asym = -1, see `tau_decay.py`); the semitauonic/baryonic
#            sources have a Dalitz-dependent tau polarization and are treated
#            unpolarized here (a full treatment is the spin-correlation item in
#            `REMAINING_WORK.md`).
TAU_SOURCES = [
    (431,  M_DS,       None, None,       '2body',  'charm',  FRAG_C[431],   BR_DS_TAUNU,        TAU_MESON_2BODY_ASYMMETRY),
    (411,  M_DPLUS,    None, None,       '2body',  'charm',  FRAG_C[411],   BR_DPLUS_TAUNU,     TAU_MESON_2BODY_ASYMMETRY),
    (521,  M_BPLUS,    None, None,       '2body',  'bottom', FRAG_B[521],   BR_BPLUS_TAUNU,     TAU_MESON_2BODY_ASYMMETRY),
    (541,  M_BC,       None, None,       '2body',  'bc',     1.0,           BR_BC_TAUNU,        TAU_MESON_2BODY_ASYMMETRY),
    (521,  M_BPLUS,    -421, M_D0,       'pseudo', 'bottom', FRAG_B[521],   BR_BP_D0_TAUNU,     0.0),
    (521,  M_BPLUS,    -423, M_DSTAR0,   'vector', 'bottom', FRAG_B[521],   BR_BP_DSTAR0_TAUNU, 0.0),
    (511,  M_B0,       -411, M_DPLUS,    'pseudo', 'bottom', FRAG_B[511],   BR_B0_DP_TAUNU,     0.0),
    (511,  M_B0,       -413, M_DSTARP,   'vector', 'bottom', FRAG_B[511],   BR_B0_DSTARP_TAUNU, 0.0),
    (531,  M_BS,       -431, M_DS,       'pseudo', 'bottom', FRAG_B[531],   BR_BS_DS_TAUNU,     0.0),
    (531,  M_BS,       -433, M_DSSTAR,   'vector', 'bottom', FRAG_B[531],   BR_BS_DSSTAR_TAUNU, 0.0),
    (5122, M_LAMBDA_B, 4122, M_LAMBDA_C, 'baryon', 'bottom', FRAG_LAMBDA_B, BR_LB_LC_TAUNU,     0.0),
]


def build_tau_pool(n_pool, rng):
    pools = {q: (sample_meson_4vectors(n_pool, q, rng=rng), get_sigma_total(q))
             for q in ('charm', 'bottom')}
    sigmas = {q: pools[q][1] for q in pools}
    shape_hnl = init_hnlcalc('Utau')  # vtau=1 so the SM-shape dBR is nonzero

    def _xsec(quark, frag, br):
        # Bc carries SIGMA_BC_PB (both charges already, no factor 2, no frag);
        # all others use 2 * sigma_quark * fragmentation * BR.
        if quark == 'bc':
            return SIGMA_BC_PB * br
        return 2.0 * sigmas[quark] * frag * br

    xsec = np.array([_xsec(quark, frag, br)
                     for (_, _, _, _, _, quark, frag, br, _) in TAU_SOURCES])
    alloc = np.maximum(1, np.round(n_pool * xsec / xsec.sum()).astype(int))

    tau4_blocks, w_blocks, a_blocks, rows = [], [], [], []
    for (parent, m_parent, dau, mdau, kind, quark, frag, br, asym), x_i, n_i in zip(
            TAU_SOURCES, xsec, alloc):
        # Bc and Lambda_b reuse the bottom FONLL (pT, y) shape with their own mass.
        shape_pool = pools['charm'][0] if quark == 'charm' else pools['bottom'][0]
        idx = rng.integers(0, n_pool, size=n_i)
        v = meson_4vec_from_kinematics(
            shape_pool['pt'][idx], shape_pool['y'][idx], shape_pool['phi'][idx], m_parent)

        if kind == '2body':
            tau4, _ = decay_2body(v['E'], v['px'], v['py'], v['pz'],
                                  m_parent, M_TAU, 0.0, rng=rng)               # d1 = tau
            label = f"{parent} -> tau nu"
        elif kind in ('pseudo', 'vector'):
            getter = (shape_hnl.get_3body_dbr_pseudoscalar if kind == 'pseudo'
                      else shape_hnl.get_3body_dbr_vector)
            dbr = getter(str(parent), str(dau), '15')
            _, tau4, _ = decay_3body_weighted_dq2dE(
                v['E'], v['px'], v['py'], v['pz'],
                m_parent, mdau, M_TAU, 0.0, dbr_expr=dbr, coupling=1.0, rng=rng)  # d2 = tau
            label = f"{parent} -> {dau} tau nu"
        elif kind == 'baryon':
            dbr = shape_hnl.get_3body_dbr_baryon(str(parent), str(dau), '15')
            _, tau4, _ = decay_3body_weighted_dq2dm122(
                v['E'], v['px'], v['py'], v['pz'],
                m_parent, mdau, M_TAU, 0.0, dbr_expr=dbr, coupling=1.0, rng=rng)  # d2 = tau
            label = f"{parent} -> {dau} tau nu"
        else:
            raise ValueError(f"unknown tau-source kind: {kind!r}")

        tau4_blocks.append(tau4)
        w_blocks.append(np.full(n_i, x_i / n_i))
        a_blocks.append(np.full(n_i, asym))
        rows.append((label, n_i, x_i))

    tau_4v = np.vstack(tau4_blocks)
    tau_w = np.concatenate(w_blocks)
    tau_asym = np.concatenate(a_blocks)

    for label, n_i, x_i in rows:
        print(f"  {label:<24s} {n_i:>6d} tau, w_sum = {x_i:.3e} pb")
    print(f"  {'Tau pool':<24s} {len(tau_w):>6d} taus, w_sum = {tau_w.sum():.3e} pb")
    return tau_4v, tau_w, tau_asym


def process_flavor(flavor, tau_4v, tau_w, tau_asym, masses, rng):
    hnl = init_hnlcalc(flavor)

    if len(tau_w) == 0:
        for m_N in masses:
            write_empty_csv(llp_csv_path(flavor, "induced_tau", m_N, base=OUTPUT_BASE))
        return

    tau_E, tau_px = tau_4v[:, 0], tau_4v[:, 1]
    tau_py, tau_pz = tau_4v[:, 2], tau_4v[:, 3]
    asym_values = np.unique(tau_asym)

    for m_N in masses:
        csv_path = llp_csv_path(flavor, "induced_tau", m_N, base=OUTPUT_BASE)
        if m_N >= M_TAU:
            write_empty_csv(csv_path)
            continue

        br_2body, br_3body, br_total = compute_tau_production_br_components(hnl, m_N)
        if br_total <= 0:
            write_empty_csv(csv_path)
            continue

        # Decay each polarization block with its own analyzing power: the 2-body
        # leptonic sources use a fixed longitudinal asymmetry; the
        # semitauonic/baryonic ones are treated unpolarized. The weight
        # (tau_w * br_total) is unaffected by the polarization; only the 2-body
        # tau->meson N angular shape changes.
        out_rows = []
        for asym in asym_values:
            sel = tau_asym == asym
            if not sel.any():
                continue
            hnl_4v, _, _ = sample_hnl_from_tau(
                tau_E[sel], tau_px[sel], tau_py[sel], tau_pz[sel], m_N,
                br_2body, br_3body, br_total, rng, asymmetry=asym,
            )
            w = tau_w[sel] * br_total
            out_rows.append(np.column_stack([
                w, hnl_4v[:, 0], hnl_4v[:, 1], hnl_4v[:, 2], hnl_4v[:, 3],
            ]))

        data = np.vstack(out_rows)
        write_csv_matrix(csv_path, data)
        print(f"    {csv_path.name}: {len(data)} events, "
              f"BR_total={br_total:.3e}, w_sum={data[:, 0].sum():.3e}")


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Generate induced-tau HNL CSVs (heavy hadron -> tau -> N)")
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

    print("Building tau pool from heavy-hadron decays ...")
    tau_4v, tau_w, tau_asym = build_tau_pool(args.n_pool, rng)

    for flavor in args.flavor:
        print(f"\n{'='*60}\nFlavor: {flavor}\n{'='*60}")
        process_flavor(flavor, tau_4v, tau_w, tau_asym, masses, rng)

    print("\nDone.")


if __name__ == "__main__":
    main()
