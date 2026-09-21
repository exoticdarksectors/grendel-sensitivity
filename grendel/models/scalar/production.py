"""BC4 production: inclusive b-hadron -> X_s S four-vectors.

Sample the bottom-hadron (pT, y, phi) pool from the FONLL spectra at 14 TeV,
rebuild the B+/B0/Bs/Lambda_b four-vectors, do the two-body decay
b-hadron -> X_s S, and write the S four-vectors weighted by
sigma_FONLL * f_frag * BR(b -> X_s S) evaluated at sin^2 theta = 1 (the
coupling is factored out and re-applied by the sensitivity scan). Inclusive
production is used because the detector reconstructs only the S decay
vertex -- the prompt X_s system is invisible, so we sum over it.

Output CSV columns are the shared four-vector format ``weight, E, px, py, pz``
(``grendel.io.vectors``); the per-event production weight is what makes the
coupling-controlled N_signal absolute.
"""
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

# Scalar mass grid for the BC4 scan: dense where the visible BR changes fastest
# (mu mu / pi pi / K K thresholds below 2 GeV), coarser up to the B -> K S
# kinematic ceiling (~4.79 GeV).
# sorted-set dedupe: np.arange float accumulation can leak a point that
# rounds onto the next segment's start (0.4999... -> 0.500).
MASS_GRID = sorted({round(x, 3) for x in (
    # Below the dimuon threshold (2 m_mu = 0.2113) the only visible channel
    # is S -> e e: no hadronic modes open below 2 m_pi and the leptonic width
    # is m_e^2-suppressed, so ctau jumps by ~3 orders of magnitude across
    # 2 m_mu (5.7 m -> 8.2 km at sin^2 theta = 1e-6).  Expect the island to
    # migrate to much larger sin^2 theta and possibly pinch off; the
    # acceptance MC uses the two-electron 2-body mode (charged fraction 1).
    list(np.arange(0.140, 0.220, 0.020)) +     # e e only (below 2 m_mu)
    list(np.arange(0.220, 0.500, 0.020)) +     # 2 m_mu -> ~K K region
    list(np.arange(0.500, 1.000, 0.025)) +
    list(np.arange(1.000, 2.000, 0.050)) +
    list(np.arange(2.000, 3.600, 0.100)) +
    list(np.arange(3.600, 4.700, 0.100)) + [4.700]
)})

# Inclusive b -> s S is a b-quark process, so every b-hadron contributes; we sum
# over the full pool: B+, B0, Bs, and the b-baryons (lumped as Lambda_b, the
# FRAG_LAMBDA_B = 0.18755 fraction the hnl set omits). The rate BR(b -> X_s S) is
# spectator-independent (model.br_B_to_Xs_S), so a baryon enters on the same
# footing as a meson -- it is only the recoil mass and lifetime that differ.
#
# The recoil mass sets the two-body S momentum and the production ceiling
# (m_S < m_parent - m_recoil). We use the *lightest strange hadron* of the right
# type as the inclusive-minimum recoil proxy: the kaon for mesons (Bs keeps the
# kaon proxy for its s-sbar recoil, a pre-existing approximation), and the Lambda
# for the b-baryon. The two choices are NOT interchangeable: at a ~5.3-5.6 GeV
# parent the S spectrum is nearly recoil-insensitive for light m_S (deep-reach
# region), but near the kinematic closure the recoil sets both the spectrum shape
# and the ceiling -- m_Lambda closes Lambda_b at 4.50 GeV, whereas a kaon proxy
# would wrongly leak it to 5.13 GeV. Using m_Lambda is therefore the physically
# correct choice, not m_K.
# PDG id -> (parent tag, m_parent, m_recoil, fragmentation frac).
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
    """S four-vectors + production weights (at sin^2 theta = 1) for mass ``m_S``.

    Returns ``(weights, E, px, py, pz)`` arrays (empty if production is closed).
    """
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
        # Two-body b-hadron -> X_s(m_recoil proxy) S(m_S); decay_2body returns (d1, d2=S).
        _, s4 = decay_2body(v["E"], v["px"], v["py"], v["pz"], m_B, m_recoil, m_S, rng=rng)
        # factor 2: b and bbar both hadronize to a b-hadron (matches hnl convention).
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
    # One shared (pT, y, phi) pool serves every mass (the boost spectrum is
    # mass-independent in the shared-shape FONLL approximation hnl also uses).
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
