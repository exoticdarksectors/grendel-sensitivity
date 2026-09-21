"""Inverse-CDF sampler from 2D FONLL dsigma/dpT/dy grids."""

import numpy as np
from . import fonll_parser
from .fonll_parser import parse_fonll_file
from ..constants import QUARK_MESON_MAP, MESON_MASSES


def _node_bin_edges(values):
    if len(values) < 2:
        raise ValueError("FONLL sampling requires at least two grid nodes")
    edges = np.empty(len(values) + 1)
    edges[0] = values[0]
    edges[-1] = values[-1]
    edges[1:-1] = 0.5 * (values[:-1] + values[1:])
    return edges


def _sample_node_intervals(edges, indices, rng):
    return rng.uniform(edges[indices], edges[indices + 1])


def _build_cdf(pt_arr, y_arr, dsigma_2d):
    n_pt = len(pt_arr)
    n_y = len(y_arr)

    pt_edges = _node_bin_edges(pt_arr)
    y_edges = _node_bin_edges(y_arr)
    pt_widths = np.diff(pt_edges)
    y_widths = np.diff(y_edges)

    prob_2d = dsigma_2d * pt_widths[:, None] * y_widths[None, :]
    prob_2d = np.maximum(prob_2d, 0.0)

    prob_flat = prob_2d.ravel()
    cdf = np.cumsum(prob_flat)
    norm = cdf[-1]
    if not np.isfinite(norm) or norm <= 0:
        raise ValueError("FONLL sampling grid has no finite positive weight")
    cdf /= norm

    return cdf, pt_edges, y_edges


def meson_4vec_from_kinematics(pt, y, phi, m):
    """Build a meson lab four-vector from sampled (pT, y, phi) and a parent mass.

    In the shared-shape approximation, one pool serves any parent species by
    re-deriving E and pz with the species mass. Returns a dict of arrays
    ``{E, px, py, pz}``.
    """
    mt = np.sqrt(np.asarray(pt) ** 2 + m ** 2)
    return {
        'E': mt * np.cosh(y),
        'px': pt * np.cos(phi),
        'py': pt * np.sin(phi),
        'pz': mt * np.sinh(y),
    }


def sample_meson_4vectors(
    n_events,
    quark,
    rng=None,
    force_species=None,
    high_pt_tilt_scale=None,
    nominal_mixture_fraction=0.5,
):
    """Sample FONLL meson kinematics, optionally with high-pT importance.

    When ``high_pt_tilt_scale`` is set, bins are drawn from a mixture of the
    nominal FONLL distribution and a distribution tilted by
    ``exp(pT / high_pt_tilt_scale)``.  ``sampling_weight`` is the exact
    nominal/proposal probability ratio for each event.  Downstream weighted
    estimators must multiply their ordinary event weight by this ratio.
    """
    if rng is None:
        rng = np.random.default_rng()

    path = fonll_parser.FONLL_FILES[quark]
    pt_arr, y_arr, dsigma_2d = parse_fonll_file(path)
    nominal_cdf, pt_edges, y_edges = _build_cdf(pt_arr, y_arr, dsigma_2d)

    n_pt = len(pt_arr)
    n_y = len(y_arr)

    if high_pt_tilt_scale is None:
        proposal_cdf = nominal_cdf
        probability_ratio = None
    else:
        scale = float(high_pt_tilt_scale)
        fraction = float(nominal_mixture_fraction)
        if not np.isfinite(scale) or scale <= 0.0:
            raise ValueError("high_pt_tilt_scale must be finite and positive")
        if not 0.0 < fraction <= 1.0:
            raise ValueError("nominal_mixture_fraction must be in (0, 1]")
        nominal = np.diff(np.concatenate([[0.0], nominal_cdf]))
        log_tilt = np.repeat(pt_arr, n_y) / scale
        log_tilt -= log_tilt.max()
        tilted = nominal * np.exp(log_tilt)
        tilted /= tilted.sum()
        proposal = fraction * nominal + (1.0 - fraction) * tilted
        proposal_cdf = np.cumsum(proposal)
        proposal_cdf[-1] = 1.0
        probability_ratio = np.divide(
            nominal,
            proposal,
            out=np.zeros_like(nominal),
            where=proposal > 0.0,
        )

    u = rng.random(n_events)
    flat_idx = np.searchsorted(proposal_cdf, u)
    flat_idx = np.clip(flat_idx, 0, n_pt * n_y - 1)

    i_pt = flat_idx // n_y
    i_y = flat_idx % n_y

    pt_sampled = _sample_node_intervals(pt_edges, i_pt, rng)
    y_sampled = _sample_node_intervals(y_edges, i_y, rng)

    phi = rng.uniform(0, 2 * np.pi, n_events)

    if force_species is not None:
        species_pdg = np.full(n_events, int(force_species), dtype=int)
        m = np.full(n_events, MESON_MASSES[int(force_species)], dtype=float)
    else:
        species_list = QUARK_MESON_MAP[quark]
        pdg_ids = np.array([s[0] for s in species_list])
        fractions = np.array([s[1] for s in species_list])
        frac_cumsum = np.cumsum(fractions)
        frac_cumsum /= frac_cumsum[-1]

        u_species = rng.random(n_events)
        species_idx = np.searchsorted(frac_cumsum, u_species)
        species_idx = np.clip(species_idx, 0, len(pdg_ids) - 1)
        species_pdg = pdg_ids[species_idx]

        m = np.array([MESON_MASSES[pdg] for pdg in species_pdg])

    mt = np.sqrt(pt_sampled**2 + m**2)
    pz = mt * np.sinh(y_sampled)
    E = mt * np.cosh(y_sampled)
    px = pt_sampled * np.cos(phi)
    py = pt_sampled * np.sin(phi)

    sampling_weight = (
        np.ones(n_events)
        if probability_ratio is None
        else probability_ratio[flat_idx]
    )
    return {
        'E': E, 'px': px, 'py': py, 'pz': pz,
        'species_pdg': species_pdg,
        'pt': pt_sampled, 'y': y_sampled, 'phi': phi,
        'sampling_weight': sampling_weight,
    }
