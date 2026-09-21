"""The analytic two-body decay backend of the dark scalar.

The first published BC4 curve decays the scalar with this engine rather
than with cached templates: per decay a mode is drawn by its branching ratio
(mu mu / pi pi / K K / s s / c c / g g / tau tau / e e), the two charged
daughters are generated back to back in the rest frame and boosted to the
lab, and the shared four-hit reconstruction and selection decide whether the
decay is reconstructed. Neutral sub-modes (pi0 pi0, K0 K0bar) give no tracks
and fail, so the visible branching fraction is folded into the outcome
exactly as the template backends fold it -- there is no separate BR_vis
multiplier.

The decay geometry is independent of sin^2 theta (every partial width scales
with it, so the mode composition does not), which is why the Monte Carlo is
built once per mass and ``scan_u2`` reweights it to every coupling.
"""
from __future__ import annotations

import numpy as np

from ...constants import CMS_ORIGIN
from ...geometry.reco_common import SIGMA_T_DEFAULT
from ...production.decay_engine.kinematics import _boost_to_lab
from ...reco.acceptance import (HIT_RESOLUTION, P_CUT, reconstruct_decays,
                                scatter_mc, selection_mask)
from . import model

# Decay mode -> (daughter mass for the 2-body proxy, charged-track fraction).
# The charged fraction encodes the neutral isospin sub-modes that give no
# tracks (pi0 pi0 = 1/3 of pi pi; K0 K0bar = 1/2 of K K). The quark/gluon
# continuum (s s, c c, g g) and 4 pi are modelled by their two LEADING charged
# hadrons (pion-mass proxy) -- a documented approximation; these channels are
# sub-dominant in the < 2 GeV region that drives the reach.
_MODE_DAUGHTER = {
    "ee":     (model.M_ELECTRON, 1.0),
    "mumu":   (model.M_MUON,     1.0),
    "tautau": (model.M_TAU,      1.0),     # tau directions proxy the leading tracks
    "pipi":   (model.M_PIPLUS,   2.0 / 3.0),
    "KK":     (model.M_KPLUS,    0.5),
    "4pi":    (model.M_PIPLUS,   1.0),
    "ss":     (model.M_PIPLUS,   1.0),
    "cc":     (model.M_PIPLUS,   1.0),
    "gg":     (model.M_PIPLUS,   1.0),
}

# Bound the daughter-to-wall ray-cast allocation. Decays are independent; the
# chunk fixes the RNG ordering and is therefore part of the published policy.
EVENT_RECO_CHUNK = 25_000


def _mode_table(m_S, width_scheme="winkler"):
    """(modes, probs, daughter_mass, charged_frac) arrays for ``m_S``."""
    br = model.branching_ratios(m_S, scheme=width_scheme)
    modes = [k for k in br if br[k] > 0.0]
    probs = np.array([br[k] for k in modes], float)
    probs = probs / probs.sum()
    m_d = np.array([_MODE_DAUGHTER[k][0] for k in modes], float)
    chf = np.array([_MODE_DAUGHTER[k][1] for k in modes], float)
    return modes, probs, m_d, chf


def build_event_mc(p4, direction, entry_d, exit_d, m_S, n_samples, rng,
                   sigma_hit=HIT_RESOLUTION, sigma_t=SIGMA_T_DEFAULT,
                   origin=CMS_ORIGIN, width_scheme="winkler",
                   reco_chunk=EVENT_RECO_CHUNK, return_mc=False):
    """Sample decay vertices + scalar decays and reconstruct, for a batch of S
    four-vectors. Returns ``(d, passed)``, each ``(n_events, n_samples)``;
    both are sin^2 theta-independent (``scan_u2`` applies the coupling). With
    ``return_mc`` additionally returns the full-length ``scatter_mc``
    selection arrays for the sequential cutflow."""
    origin = np.asarray(origin, float)
    p4 = np.asarray(p4, float)
    direction = np.asarray(direction, float)
    n_ev = len(entry_d)

    d = rng.uniform(np.asarray(entry_d)[:, None], np.asarray(exit_d)[:, None],
                    size=(n_ev, n_samples))
    M = n_ev * n_samples
    vtx = origin[None, :] + d.reshape(M, 1) * np.repeat(direction, n_samples, axis=0)
    parent = np.repeat(p4, n_samples, axis=0)            # (M, 4) S four-vectors

    # Sample a decay mode per decay, then its charged/neutral sub-mode.
    modes, probs, m_d_tab, chf_tab = _mode_table(m_S, width_scheme=width_scheme)
    mi = rng.choice(len(modes), size=M, p=probs)
    m_d = m_d_tab[mi]
    charged = rng.random(M) < chf_tab[mi]
    pstar = np.sqrt(np.clip(0.25 * m_S ** 2 - m_d ** 2, 0.0, None))
    valid = charged & (pstar > 0)

    passed = np.zeros(M, dtype=bool)
    mc_full = None
    if valid.any():
        idx = np.where(valid)[0]
        # Isotropic back-to-back daughters in the S rest frame.
        cos_t = rng.uniform(-1.0, 1.0, len(idx))
        sin_t = np.sqrt(np.clip(1.0 - cos_t ** 2, 0.0, None))
        phi = rng.uniform(0.0, 2 * np.pi, len(idx))
        nhat = np.column_stack([sin_t * np.cos(phi), sin_t * np.sin(phi), cos_t])
        ps = pstar[idx][:, None]
        E_d = np.full(len(idx), m_S / 2.0)
        d1_rest = np.column_stack([E_d, ps * nhat])
        d2_rest = np.column_stack([E_d, -ps * nhat])
        pe, ppx, ppy, ppz = (parent[idx, 0], parent[idx, 1],
                             parent[idx, 2], parent[idx, 3])
        d1 = _boost_to_lab(d1_rest, pe, ppx, ppy, ppz)
        d2 = _boost_to_lab(d2_rest, pe, ppx, ppy, ppz)
        p1 = np.linalg.norm(d1[:, 1:], axis=1)
        p2 = np.linalg.norm(d2[:, 1:], axis=1)
        dir1 = d1[:, 1:] / p1[:, None]
        dir2 = d2[:, 1:] / p2[:, None]
        p_soft = np.minimum(p1, p2)
        # Daughter speeds beta = |p|/E feed the timing model (soft pions and
        # kaons are genuinely slow; beta = 1 would overestimate the timing-cut
        # acceptance). Guard against E < |p| roundoff after the boost.
        b1 = p1 / np.maximum(d1[:, 0], p1)
        b2 = p2 / np.maximum(d2[:, 0], p2)
        # Both daughters must clear the track momentum floor (mirrors best-two).
        ok = (p_soft > P_CUT)
        if ok.any():
            sub = idx[ok]
            ok_idx = np.where(ok)[0]
            if return_mc:
                mc_full = scatter_mc(None, np.zeros(M, dtype=bool), M)
                mc_full['visible'][sub] = True
            for start in range(0, len(sub), reco_chunk):
                stop = min(start + reco_chunk, len(sub))
                take = ok_idx[start:stop]
                mc = reconstruct_decays(
                    vtx[sub[start:stop]], dir1[take], dir2[take], p_soft[take],
                    sigma_hit, sigma_t, rng, beta1=b1[take], beta2=b2[take])
                passed[sub[start:stop]] = selection_mask(mc)
                if return_mc:
                    rows_idx = sub[start:stop]
                    for key in ('sep', 'sep_outer', 'open_angle', 'dca',
                                'collin', 'pointing', 'vtx_in', 'on_tracker',
                                'timing_chi2', 'p_soft'):
                        mc_full[key][rows_idx] = mc[key]

    if return_mc:
        if mc_full is None:
            mc_full = scatter_mc(None, np.zeros(M, dtype=bool), M)
        return d, passed.reshape(n_ev, n_samples), mc_full
    return d, passed.reshape(n_ev, n_samples)


class AnalyticBackend:
    """``DecayBackend`` for the two-body proxy above."""

    def __init__(self, m_S, width_scheme="winkler"):
        self.m_S = float(m_S)
        self.width_scheme = width_scheme
        self.name = f"analytic-2body:{width_scheme}"
        self.ctau_ref = model.ctau_sin2theta1(self.m_S, scheme=width_scheme)

    def build(self, p4, direction, entry_d, exit_d, n_samples, rng, return_mc=False):
        out = build_event_mc(p4, direction, entry_d, exit_d, self.m_S, n_samples, rng,
                             width_scheme=self.width_scheme, return_mc=return_mc)
        if return_mc:
            d, passed, mc = out
            return d, passed, None, mc
        d, passed = out
        return d, passed, None, None

    def sample_weights(self, template_index):
        return None
