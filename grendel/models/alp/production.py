#!/usr/bin/env python3
"""BC10 production: B -> K^(*) a four-vectors at the reference coupling.

Sample the bottom-hadron (pT, y, phi) pool from the FONLL spectra at 14 TeV,
rebuild the B+ and B0 four-vectors, decay them two-body to a kaon of the
kaon tower (K, K0*, K*, K1, K2*) plus the ALP, and write the ALP four-vectors
weighted by sigma_FONLL * f_frag * BR(B -> K a) at 1/f = 1/f_ref. The
coupling is factored out and re-applied by the sensitivity scan. Output
columns are the shared four-vector format ``weight, E, px, py, pz``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from ...io.paths import ModelPaths
from ...io.vectors import format_mass_for_filename, write_empty_csv, write_llp_csv
from ...production.constants import FRAG_B
from ...production.decay_engine.kinematics import decay_2body
from ...production.fonll.fonll_parser import get_sigma_total
from ...production.fonll.meson_sampler import (meson_4vec_from_kinematics,
                                                sample_meson_4vectors)
from . import model
from .mass_grid import ALP_MASS_GRID


def alp_csv_path(mass, base):
    path = Path(base) / f"mA_{format_mass_for_filename(mass)}.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _atomic_write(path, weights=None, energy=None, px=None, py=None, pz=None):
    """Write one mass output atomically so ``--resume`` never accepts a partial CSV."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    if weights is None:
        write_empty_csv(temporary)
    else:
        write_llp_csv(temporary, weights, energy, px, py, pz)
    temporary.replace(path)


def _importance_weighted_rate(normalization, pool, indices):
    """Apply the nominal/proposal FONLL density ratio to selected parents."""
    ratios = pool.get("sampling_weight")
    if ratios is None:
        return np.full(len(indices), normalization)
    return normalization * np.asarray(ratios, dtype=float)[indices]


def _partition_parent_indices(n_pool, n_channels, rng):
    """Assign disjoint, random parent samples to each production channel."""
    if n_pool < n_channels:
        raise ValueError("n_pool must be at least the number of channels")
    n_each = n_pool // n_channels
    usable = n_each * n_channels
    permutation = rng.permutation(n_pool)[:usable]
    return permutation.reshape(n_channels, n_each)


def generate(
    masses,
    n_pool,
    out_dir,
    seed=42,
    resume=False,
    cbs_amplitude_scale=1.0,
    high_pt_tilt_scale=None,
    nominal_mixture_fraction=0.5,
):
    if cbs_amplitude_scale <= 0.0:
        raise ValueError("cbs_amplitude_scale must be positive")
    rng = np.random.default_rng(seed)
    print(f"Sampling {n_pool} FONLL bottom (pt,y,phi) shape events...", flush=True)
    pool = sample_meson_4vectors(
        n_pool,
        "bottom",
        rng=rng,
        high_pt_tilt_scale=high_pt_tilt_scale,
        nominal_mixture_fraction=nominal_mixture_fraction,
    )
    sigma_b = get_sigma_total("bottom")
    print(f"  sigma_FONLL(bottom) = {sigma_b:.4e} pb", flush=True)

    channels = [(label, pdg, kaon)
                for label, pdg in model.PRODUCTION_PARENTS
                for kaon in model.KAON_TOWER]

    for m_a in masses:
        resonance = model.excluded_light_meson_resonance(m_a)
        if resonance is not None:
            print(
                f"  m_a={m_a:.3f}: skip unsupported {resonance} resonance",
                flush=True,
            )
            continue
        path = alp_csv_path(m_a, out_dir)
        already_done = resume and path.exists()
        all_w, all_E, all_px, all_py, all_pz = [], [], [], [], []
        parent_indices = _partition_parent_indices(
            n_pool, len(channels), rng,
        )
        for channel_index, (label, pdg, kaon) in enumerate(channels):
            m_B = model.M_BPLUS if label == "B+" else model.M_B0
            m_K = model.kaon_mass(kaon, label)
            if m_a >= m_B - m_K:
                continue
            br = model.br_B_to_Ka(
                m_a, model.INV_F_REF, parent=label, kaon=kaon
            ) * cbs_amplitude_scale**2
            if br <= 0.0:
                continue
            frag = FRAG_B[pdg]
            idx = parent_indices[channel_index]
            n_each = len(idx)
            v = meson_4vec_from_kinematics(
                pool["pt"][idx], pool["y"][idx], pool["phi"][idx], m_B)
            # two-body B -> K_i(m_K) + a(m_a); decay_2body returns (d1=K, d2=a)
            _, a4 = decay_2body(v["E"], v["px"], v["py"], v["pz"],
                                m_B, m_K, m_a, rng=rng)
            w = _importance_weighted_rate(
                2.0 * sigma_b * frag * br / n_each,
                pool,
                idx,
            )
            if not already_done:
                all_w.append(w)
                all_E.append(a4[:, 0]); all_px.append(a4[:, 1])
                all_py.append(a4[:, 2]); all_pz.append(a4[:, 3])

        if already_done:
            print(f"  m_a={m_a:.3f}: already checkpointed -> {path.name}", flush=True)
            continue
        if all_w:
            weights = np.concatenate(all_w)
            _atomic_write(
                path, weights, np.concatenate(all_E), np.concatenate(all_px),
                np.concatenate(all_py), np.concatenate(all_pz),
            )
            print(f"  m_a={m_a:.3f}: {len(weights)} a, "
                  f"sigma_ref={weights.sum():.3e} pb -> {path.name}", flush=True)
        else:
            _atomic_write(path)
            print(f"  m_a={m_a:.3f}: 0 (above B->K a threshold)", flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description="BC10 a four-vector production")
    ap.add_argument("--mass", type=float, nargs="+", default=None)
    ap.add_argument("--n-pool", type=int, default=200_000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument(
        "--resume", action="store_true",
        help="keep atomically completed mass CSVs while replaying the RNG stream",
    )
    ap.add_argument(
        "--cbs-amplitude-scale", type=float, default=1.0,
        help="multiply C_bs by this factor (production rates scale as its square)",
    )
    ap.add_argument(
        "--high-pt-tilt-scale",
        type=float,
        default=None,
        help=(
            "importance-sample a nominal/high-pT FONLL mixture tilted by "
            "exp(pT/scale); output weights include the exact p/q correction"
        ),
    )
    ap.add_argument(
        "--nominal-mixture-fraction",
        type=float,
        default=0.5,
        help="nominal fraction of the optional high-pT proposal mixture",
    )
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="four-vector directory (default: $GRENDEL_BC10_VECTORS_DIR)")
    args = ap.parse_args(argv)
    masses = args.mass if args.mass else ALP_MASS_GRID
    out_dir = args.out_dir or ModelPaths.resolve("bc10").vectors
    print(f"BC10 production: {len(masses)} masses, n_pool={args.n_pool}")
    print(f"  output -> {out_dir}")
    generate(
        masses, args.n_pool, out_dir, args.seed,
        resume=args.resume,
        cbs_amplitude_scale=args.cbs_amplitude_scale,
        high_pt_tilt_scale=args.high_pt_tilt_scale,
        nominal_mixture_fraction=args.nominal_mixture_fraction,
    )
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
