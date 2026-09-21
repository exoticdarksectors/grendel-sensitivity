#!/usr/bin/env python3
"""Charged-kaon -> HNL production CSVs with a Pythia SoftQCD spectrum + transport.

Default model (Pythia spectrum + transport survival weight):
  * spectrum   -- Pythia 8.315 SoftQCD:inelastic (pT, y) histogram, committed in
                  ``production/data/kaon_softqcd_spectrum.npz`` (built by
                  ``kaon_softqcd.cc``, rebuilt from tracked sources via
                  ``make_kaon_spectrum.py``); normalization sigma_inel * <n_K+->.
  * transport  -- a charged kaon (ctau = 3.712 m) is absorbed in dense material
                  unless it decays first; each kaon carries a survival weight
                  P(decay within d_esc) = 1 - exp(-d_esc / (beta*gamma * ctau)).
                  The displaced decay origin shifts the HNL start by <= d_esc, which
                  make_kaon_spectrum's sibling transport_control.py measures to change
                  the acceptance by only a few % on the sensitivity plateau, so the
                  HNL is cast from IP. ``d_esc`` (KAON_D_ESC = 1.5 m) is varied over
                  KAON_D_ESC_RANGE = [1, 3] m in the published BC6/BC7 transport band;
                  it remains a proxy pending a CMS material map.

Legacy behaviour (``--spectrum tsallis`` and/or ``--no-transport``) reproduces
the old parametric prompt-at-IP stub for A/B comparison.
"""

import argparse
import random
import sys
from pathlib import Path

import numpy as np


from grendel.models.hnl.mass_grid import MASS_GRID
from grendel.production.constants import (
    M_KAON, LEPTON_MASSES, FLAVOR_TO_LEPTON_PDG,
    SIGMA_KAON_PB, SIGMA_KAON_PB_TSALLIS, KAON_D_ESC,
    KAON_TSALLIS_T, KAON_TSALLIS_N, KAON_PT_MAX, KAON_RAPIDITY_SIGMA, KAON_E_MAX,
)
from grendel.models.hnl.production.generate_meson_csvs import (
    compute_production_br_components, _sample_hnl_from_mesons,
)
from grendel.models.hnl.production.hnlcalc import init_hnlcalc
from grendel.models.hnl.production.paths import llp_csv_path, write_empty_csv, write_llp_csv
from grendel.models.hnl.production.paths import LLP_VECTORS_DIR

OUTPUT_BASE = LLP_VECTORS_DIR
N_POOL = 100_000
KAON_PDG = 321
CTAU_KAON = 3.712                       # m, c*tau(K+-)
_SPECTRUM_NPZ = Path(__file__).resolve().parent / "data" / "kaon_softqcd_spectrum.npz"


def _kaon_4vectors_from_ptphi(pt, y, phi):
    mt = np.sqrt(pt * pt + M_KAON * M_KAON)
    return {
        "E": mt * np.cosh(y), "px": pt * np.cos(phi),
        "py": pt * np.sin(phi), "pz": mt * np.sinh(y),
        "species_pdg": np.full(len(pt), KAON_PDG, dtype=int),
        "pt": pt, "y": y, "phi": phi,
    }


def sample_kaon_4vectors_pythia(n_pool, rng, spectrum_path=_SPECTRUM_NPZ):
    """Sample (pT, y) from the committed Pythia SoftQCD histogram (inverse CDF on
    the flattened 2D histogram), then rebuild the K+- four-vectors."""
    data = np.load(spectrum_path)
    hist = data["hist"]
    pt_edges, y_edges = data["pt_edges"], data["y_edges"]
    flat = hist.ravel().astype(float)
    cdf = np.cumsum(flat)
    cdf /= cdf[-1]
    idx = np.searchsorted(cdf, rng.random(n_pool))
    idx = np.clip(idx, 0, len(flat) - 1)
    i_pt, i_y = np.unravel_index(idx, hist.shape)
    pt = rng.uniform(pt_edges[i_pt], pt_edges[i_pt + 1])
    y = rng.uniform(y_edges[i_y], y_edges[i_y + 1])
    phi = rng.uniform(0.0, 2.0 * np.pi, n_pool)
    return _kaon_4vectors_from_ptphi(pt, y, phi)


def sample_kaon_4vectors_tsallis(n_pool, rng):
    """Legacy parametric spectrum: Tsallis pT x Gaussian rapidity."""
    m_k = M_KAON
    n_grid = 4000
    pt_edges = np.linspace(0.0, KAON_PT_MAX, n_grid + 1)
    pt_centers = 0.5 * (pt_edges[:-1] + pt_edges[1:])
    mt = np.sqrt(pt_centers**2 + m_k**2)
    pdf = pt_centers * (1.0 + (mt - m_k) / (KAON_TSALLIS_N * KAON_TSALLIS_T)) ** (-KAON_TSALLIS_N)
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]
    idx = np.clip(np.searchsorted(cdf, rng.random(n_pool)), 0, n_grid - 1)
    pt = rng.uniform(pt_edges[idx], pt_edges[idx + 1])
    y = rng.normal(0.0, KAON_RAPIDITY_SIGMA, n_pool)
    mt = np.sqrt(pt * pt + m_k * m_k)
    bad = mt * np.cosh(y) >= KAON_E_MAX
    while bad.any():
        y[bad] = rng.normal(0.0, KAON_RAPIDITY_SIGMA, int(bad.sum()))
        bad = mt * np.cosh(y) >= KAON_E_MAX
    phi = rng.uniform(0.0, 2.0 * np.pi, n_pool)
    return _kaon_4vectors_from_ptphi(pt, y, phi)


def sample_kaon_4vectors(n_pool, rng, spectrum="pythia"):
    if spectrum == "pythia":
        return sample_kaon_4vectors_pythia(n_pool, rng)
    if spectrum == "tsallis":
        return sample_kaon_4vectors_tsallis(n_pool, rng)
    raise ValueError(f"unknown kaon spectrum {spectrum!r}")


def kaon_survival_weight(pool, d_esc):
    """Per-kaon transport survival P(decay within d_esc before absorption) =
    1 - exp(-d_esc / (beta*gamma * ctau_K)). d_esc=None -> prompt (all 1.0)."""
    n = len(pool["E"])
    if d_esc is None:
        return np.ones(n)
    p = np.sqrt(pool["px"]**2 + pool["py"]**2 + pool["pz"]**2)
    beta_gamma = p / M_KAON
    decay_len = beta_gamma * CTAU_KAON
    return 1.0 - np.exp(-d_esc / decay_len)


def process_flavor(flavor, pool, masses, rng,
                   sigma_kaon_pb=SIGMA_KAON_PB, d_esc=KAON_D_ESC):
    hnl = init_hnlcalc(flavor)
    lepton_pdg = FLAVOR_TO_LEPTON_PDG[flavor]
    m_lepton = LEPTON_MASSES[flavor]

    n_pool = len(pool["E"])
    m_parent = M_KAON
    survival = kaon_survival_weight(pool, d_esc)   # per-kaon transport weight

    for m_N in masses:
        csv_path = llp_csv_path(flavor, "Kmeson", m_N, base=OUTPUT_BASE)
        if m_N >= m_parent - m_lepton:
            write_empty_csv(csv_path)
            continue

        _, br_3body_channels, br = compute_production_br_components(
            hnl, KAON_PDG, lepton_pdg, m_N)
        if br <= 0:
            write_empty_csv(csv_path)
            continue

        hnl_4v = _sample_hnl_from_mesons(
            pool["E"], pool["px"], pool["py"], pool["pz"],
            m_parent, m_lepton, m_N, br_3body_channels, br, hnl, rng)
        # per-kaon weight: flux * BR / n_pool, scaled by the transport survival
        weights = (sigma_kaon_pb * br / n_pool) * survival
        write_llp_csv(csv_path, weights, hnl_4v[:, 0], hnl_4v[:, 1],
                      hnl_4v[:, 2], hnl_4v[:, 3])
        print(f"    {csv_path.name}: {n_pool} events, BR={br:.3e}, "
              f"w_sum={weights.sum():.3e} (survival<w>={survival.mean():.3f})")


def main():
    parser = argparse.ArgumentParser(description="Generate kaon -> HNL CSVs (K+ -> l N)")
    parser.add_argument("--flavor", choices=["Ue", "Umu", "Utau"], nargs="+",
                        default=["Ue", "Umu", "Utau"])
    parser.add_argument("--n-pool", type=int, default=N_POOL)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--masses", type=float, nargs="+", default=None)
    parser.add_argument("--spectrum", choices=["pythia", "tsallis"], default="pythia",
                        help="kaon (pT,y) spectrum (default: committed Pythia SoftQCD)")
    parser.add_argument("--d-esc", type=float, default=KAON_D_ESC,
                        help="transport escape distance in m (default KAON_D_ESC)")
    parser.add_argument("--no-transport", action="store_true",
                        help="legacy prompt-at-IP: no kaon material-survival weight")
    args = parser.parse_args()

    masses = args.masses if args.masses else MASS_GRID
    d_esc = None if args.no_transport else args.d_esc
    sigma_kaon_pb = SIGMA_KAON_PB if args.spectrum == "pythia" else SIGMA_KAON_PB_TSALLIS

    print(f"Kaon pool: spectrum={args.spectrum}, sigma_kaon={sigma_kaon_pb:.3e} pb, "
          f"transport d_esc={d_esc} m ({args.n_pool} events)")
    # Re-seed per flavour. This RNG must not chain across the loop: it used to
    # be seeded once above, so a flavour's sample depended on which other
    # flavours were requested and in what order, and `--flavor Ue Umu Utau` did
    # not reproduce three single-flavour invocations. Re-sampling the pool
    # inside the loop is deliberate -- it restores the RNG to the same
    # post-pool state a standalone invocation sees, so both call styles agree
    # bit-for-bit.
    for flavor in args.flavor:
        print(f"\n{'='*60}\nFlavor: {flavor}\n{'='*60}")
        rng = np.random.default_rng(args.seed)
        random.seed(args.seed)
        pool = sample_kaon_4vectors(args.n_pool, rng, spectrum=args.spectrum)
        process_flavor(flavor, pool, masses, rng, sigma_kaon_pb, d_esc)

    print("\nDone.")


if __name__ == "__main__":
    main()
