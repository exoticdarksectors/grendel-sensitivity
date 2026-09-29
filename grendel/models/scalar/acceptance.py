"""The analytic two-body decay backend of the dark scalar."""
from __future__ import annotations

import numpy as np

from ...constants import CMS_ORIGIN
from ...geometry.reco_common import SIGMA_T_DEFAULT
from ...production.decay_engine.kinematics import _boost_to_lab
from ...reco.acceptance import (HIT_RESOLUTION, P_CUT, reconstruct_decays,
                                scatter_mc, selection_mask)
from . import model

_MODE_DAUGHTER = {
    "ee":     (model.M_ELECTRON, 1.0),
    "mumu":   (model.M_MUON,     1.0),
    "tautau": (model.M_TAU,      1.0),
    "pipi":   (model.M_PIPLUS,   2.0 / 3.0),
    "KK":     (model.M_KPLUS,    0.5),
    "4pi":    (model.M_PIPLUS,   1.0),
    "ss":     (model.M_PIPLUS,   1.0),
    "cc":     (model.M_PIPLUS,   1.0),
    "gg":     (model.M_PIPLUS,   1.0),
}

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
    """Sample decay vertices + scalar decays and reconstruct, for a batch of S four-vectors."""
    origin = np.asarray(origin, float)
    p4 = np.asarray(p4, float)
    direction = np.asarray(direction, float)
    n_ev = len(entry_d)

    d = rng.uniform(np.asarray(entry_d)[:, None], np.asarray(exit_d)[:, None],
                    size=(n_ev, n_samples))
    M = n_ev * n_samples
    vtx = origin[None, :] + d.reshape(M, 1) * np.repeat(direction, n_samples, axis=0)
    parent = np.repeat(p4, n_samples, axis=0)

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
        b1 = p1 / np.maximum(d1[:, 0], p1)
        b2 = p2 / np.maximum(d2[:, 0], p2)
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
