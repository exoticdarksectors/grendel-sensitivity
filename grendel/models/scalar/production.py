"""BC4 production: inclusive b-hadron -> X_s S four-vectors."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ...production.constants import FRAG_B, FRAG_LAMBDA_B
from ...production.fonll.fonll_parser import get_sigma_total
from ...production.fonll.meson_sampler import (
    sample_meson_4vectors, meson_4vec_from_kinematics)
from ...production.decay_engine.kinematics import decay_2body
from ...io.vectors import format_mass_for_filename
from . import model

MASS_GRID = sorted({round(x, 3) for x in (
    list(np.arange(0.140, 0.220, 0.020)) +
    list(np.arange(0.220, 0.500, 0.020)) +
    list(np.arange(0.500, 1.000, 0.025)) +
    list(np.arange(1.000, 2.000, 0.050)) +
    list(np.arange(2.000, 3.600, 0.100)) +
    list(np.arange(3.600, 4.700, 0.100)) + [4.700]
)})

B_SPECIES = {
    521:  ("B+",       model.M_BPLUS,    model.M_KPLUS,  FRAG_B[521]),
    511:  ("B0",       model.M_B0,       model.M_K0,     FRAG_B[511]),
    531:  ("Bs",       model.M_BS,       model.M_KPLUS,  FRAG_B[531]),
    5122: ("Lambda_b", model.M_LAMBDA_B, model.M_LAMBDA, FRAG_LAMBDA_B),
}

N_POOL_DEFAULT = 200_000


def generate_scalar_4vectors(
    m_S,
    n_pool,
    rng,
    sigma_bottom=None,
    pool=None,
    high_pt_tilt_scale=None,
    nominal_mixture_fraction=0.5,
):
    """S four-vectors + production weights (at sin^2 theta = 1) for mass ``m_S``."""
    if m_S >= model.M_S_MAX_BTOK:
        return (np.empty(0),) * 5
    if sigma_bottom is None:
        sigma_bottom = get_sigma_total("bottom")
    if pool is None:
        pool = sample_meson_4vectors(
            n_pool,
            "bottom",
            rng=rng,
            high_pt_tilt_scale=high_pt_tilt_scale,
            nominal_mixture_fraction=nominal_mixture_fraction,
        )
    n_each = len(pool["pt"])
    sampling_weight = pool.get("sampling_weight", np.ones(n_each))

    weights, E, px, py, pz = [], [], [], [], []
    for pdg, (parent, m_B, m_recoil, frag) in B_SPECIES.items():
        if m_S >= m_B - m_recoil:
            continue
        br = float(model.br_B_to_Xs_S(m_S, parent=parent, sin2theta=1.0))
        if br <= 0:
            continue
        v = meson_4vec_from_kinematics(pool["pt"], pool["y"], pool["phi"], m_B)
        _, s4 = decay_2body(v["E"], v["px"], v["py"], v["pz"], m_B, m_recoil, m_S, rng=rng)
        w = 2.0 * sigma_bottom * frag * br / n_each
        weights.append(w * sampling_weight)
        E.append(s4[:, 0]); px.append(s4[:, 1]); py.append(s4[:, 2]); pz.append(s4[:, 3])

    if not weights:
        return (np.empty(0),) * 5
    return (np.concatenate(weights), np.concatenate(E),
            np.concatenate(px), np.concatenate(py), np.concatenate(pz))


def _mass_label(m_S):
    return format_mass_for_filename(m_S)


def write_scalar_csv(m_S, out_dir, n_pool, rng, sigma_bottom=None, pool=None):
    """Generate and write the S four-vector CSV for ``m_S`` (hnl format)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"mS_{_mass_label(m_S)}.csv"
    w, E, px, py, pz = generate_scalar_4vectors(
        m_S, n_pool, rng, sigma_bottom=sigma_bottom, pool=pool)
    if len(w) == 0:
        path.write_text("")
    else:
        np.savetxt(path, np.column_stack([w, E, px, py, pz]),
                   delimiter=",", fmt="%.8e")
    return path, len(w)


def main(argv=None):
    import argparse

    p = argparse.ArgumentParser(
        description="Generate inclusive BC4 b -> X_s S events with a K-recoil proxy"
    )
    p.add_argument("--out-dir", default=None,
                   help="four-vector directory (default: $GRENDEL_BC4_VECTORS_DIR)")
    p.add_argument("--n-pool", type=int, default=N_POOL_DEFAULT)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--masses", type=float, nargs="+", default=None)
    p.add_argument(
        "--high-pt-tilt-scale",
        type=float,
        default=None,
        help=(
            "importance-sample a nominal/high-pT FONLL mixture tilted by "
            "exp(pT/scale); output weights include the exact p/q correction"
        ),
    )
    p.add_argument(
        "--nominal-mixture-fraction",
        type=float,
        default=0.5,
        help="nominal fraction of the optional high-pT proposal mixture",
    )
    args = p.parse_args(argv)

    from ...io.paths import ModelPaths
    out_dir = args.out_dir or ModelPaths.resolve("bc4").vectors
    rng = np.random.default_rng(args.seed)
    masses = args.masses if args.masses else MASS_GRID
    sigma_bottom = get_sigma_total("bottom")
    print(f"sigma_FONLL(bottom) = {sigma_bottom:.3e} pb;  {len(masses)} masses")
    pool = sample_meson_4vectors(
        args.n_pool,
        "bottom",
        rng=rng,
        high_pt_tilt_scale=args.high_pt_tilt_scale,
        nominal_mixture_fraction=args.nominal_mixture_fraction,
    )
    for m_S in masses:
        path, n = write_scalar_csv(m_S, out_dir, args.n_pool, rng,
                                   sigma_bottom=sigma_bottom, pool=pool)
        print(f"  m_S={m_S:.3f}: {n} events -> {path.name}")


if __name__ == "__main__":
    main()
