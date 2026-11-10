"""GRENDEL signal acceptance from FairShip HNL decays (Stage 2).

Replaces the earlier idealized analytic two-body acceptance (now removed) with
the single-source PR #13 reconstruction (``reco_common``): each HNL decay's
charged daughters (from the FairShip rest-frame templates, boosted to the lab)
are turned into real wall hits and run through the same bounded 4-hit
reconstruction + selection the cosmic-decay background uses.

Pipeline per HNL four-vector that hits the fiducial air volume:
  1. sample a decay vertex along the flight path (uniform; decay-density weight
     folded in by the lifetime reweighting, as in the higgs MC);
  2. draw a FairShip rest-frame decay template and boost it to the lab;
  3. keep charged stable daughters with |p| > P_CUT, take the two highest-momentum
     (best-two-track); fewer than two -> unreconstructable;
  4. ray-cast both to the wall, require both on a tracker surface, build the
     inner/outer hits (``wall_inner_outer``), reconstruct (``reconstruct_3d``),
     and time them (``timing_chi2_4hit``);
  5. apply ``selection_mask`` (gate / pointing / collinearity / timing).

Geometry + reconstruction are the *shared* ``higgs/grendel_geometry`` /
``higgs/reco_common`` single source on exoticdarksectors/llpatcolliders main
(PR #13) -- imported directly, not copied. Selection constants are mirrored
from higgs/decayProbPerEvent_2body.py so signal and the cosmic background share
one definition.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

_HNL_ROOT = Path(__file__).resolve().parent.parent
# Single source for detector geometry + 4-hit reconstruction: import the SAME
# grendel_geometry / reco_common the higgs analysis uses (PR #13 on main), so
# signal and the cosmic background share one reconstruction definition.
_RECO_DIR = _HNL_ROOT.parent / "higgs"
if str(_RECO_DIR) not in sys.path:
    sys.path.insert(0, str(_RECO_DIR))

from grendel_geometry import (  # noqa: E402
    mesh_fiducial, points_on_tracker, DETECTOR_THICKNESS,
)
import reco_common as _rc  # noqa: E402
from reco_common import SIGMA_T_DEFAULT, CHI2_TIMING_MAX  # noqa: E402

# --- Selection constants (mirror higgs/decayProbPerEvent_2body.py on main) ---
P_CUT = float(os.environ.get("HNL_P_CUT", "0.600"))  # GeV/c — min charged-track momentum (env override for cut studies)
SEP_MIN = 0.01                # m
SEP_MAX = 10.0                # m
DCA_CUT = 0.1                 # m
THETA_PARALLEL = 0.100        # rad
SEP_OUT_MAX_PARALLEL = 0.30   # m
COLLIN_FRAC = 0.40            # collinearity > COLLIN_FRAC * L
SEP_OUT_GATE = 0.16           # m
POINT_TIGHT_SEP_IN = 0.050    # rad
SEP_IN_POINT_GATE = 0.10      # m
POINT_GLOBAL = 0.8            # rad
HIT_RESOLUTION = 0.003        # m (per-layer)


def _first_forward_hit(mesh, origins, dirs):
    """Nearest forward mesh intersection per ray, NaN on miss (from main)."""
    out = np.full((len(origins), 3), np.nan)
    locs, ray_idx, _ = mesh.ray.intersects_location(
        ray_origins=origins, ray_directions=dirs)
    if len(locs) == 0:
        return out
    signed = np.einsum('ij,ij->i', locs - origins[ray_idx], dirs[ray_idx])
    fwd = signed > 1e-6
    locs, ray_idx, signed = locs[fwd], ray_idx[fwd], signed[fwd]
    if len(locs) == 0:
        return out
    order = np.lexsort((signed, ray_idx))
    locs, ray_idx = locs[order], ray_idx[order]
    _, first = np.unique(ray_idx, return_index=True)
    out[ray_idx[first]] = locs[first]
    return out


def boost_rest_to_lab(parent_p4, rest_px, rest_py, rest_pz, rest_E):
    """General Lorentz boost of rest-frame daughter 3-momenta into the lab,
    along the parent 3-momentum (replicates fairship_decay.boost_decay_to_lab
    on flat arrays). parent_p4 = (E, px, py, pz). Returns (px, py, pz, E) lab."""
    E0 = float(parent_p4[0])
    beta = np.asarray(parent_p4[1:], float) / E0
    beta2 = float(beta @ beta)
    rest_p = np.column_stack([rest_px, rest_py, rest_pz]).astype(float)
    rest_E = np.asarray(rest_E, float)
    if beta2 <= 0.0:
        return rest_p[:, 0], rest_p[:, 1], rest_p[:, 2], rest_E
    if beta2 >= 1.0:
        beta2 = np.nextafter(1.0, 0.0)
    gamma = 1.0 / np.sqrt(1.0 - beta2)
    bp = rest_p @ beta
    g2 = (gamma - 1.0) / beta2
    lab_p = rest_p + (g2 * bp + gamma * rest_E)[:, None] * beta[None, :]
    lab_E = gamma * (rest_E + bp)
    return lab_p[:, 0], lab_p[:, 1], lab_p[:, 2], lab_E


def best_two_directions(px, py, pz, charge, stable, p_cut=P_CUT):
    """From one decay's lab daughters, pick the two highest-momentum charged
    stable tracks above p_cut. Returns (dir1, dir2, p_soft) or None if <2."""
    p = np.sqrt(px**2 + py**2 + pz**2)
    sel = (np.abs(charge) > 0.5) & stable.astype(bool) & (p > p_cut)
    if int(sel.sum()) < 2:
        return None
    idx = np.where(sel)[0]
    order = idx[np.argsort(p[idx])[::-1]]   # descending momentum
    i1, i2 = order[0], order[1]
    d1 = np.array([px[i1], py[i1], pz[i1]]) / p[i1]
    d2 = np.array([px[i2], py[i2], pz[i2]]) / p[i2]
    return d1, d2, float(min(p[i1], p[i2]))


def reconstruct_decays(decay_pos, dir1, dir2, p_soft, sigma_hit, sigma_t, rng,
                       mesh=mesh_fiducial):
    """Build the selection ``mc`` dict for N decays from their vertex + the two
    daughter directions (the shared geometry->reco core, mirroring higgs
    sample_separations lines 619-736). All inputs (N,3)/(N,)."""
    decay_pos = np.asarray(decay_pos, float)
    dir1 = np.asarray(dir1, float)
    dir2 = np.asarray(dir2, float)
    n = len(decay_pos)

    origins = np.concatenate([decay_pos, decay_pos])
    dirs = np.concatenate([dir1, dir2])
    exit_pts = _first_forward_hit(mesh, origins, dirs)
    on_trk = np.zeros(2 * n, dtype=bool)
    valid = ~np.isnan(exit_pts[:, 0])
    if valid.any():
        on_trk[valid] = points_on_tracker(exit_pts[valid])
    on_tracker = on_trk[:n] & on_trk[n:]

    sep = np.full(n, np.nan); sep_outer = np.full(n, np.nan)
    open_angle = np.full(n, np.nan); dca = np.full(n, np.nan)
    collin = np.full(n, np.nan); pointing = np.full(n, np.nan)
    vtx_in = np.zeros(n, dtype=bool)
    timing = np.full(n, np.inf)

    in1, out1 = _rc.wall_inner_outer(exit_pts[:n], dir1, DETECTOR_THICKNESS)
    in2, out2 = _rc.wall_inner_outer(exit_pts[n:], dir2, DETECTOR_THICKNESS)
    fin = (np.all(np.isfinite(in1), 1) & np.all(np.isfinite(out1), 1)
           & np.all(np.isfinite(in2), 1) & np.all(np.isfinite(out2), 1))
    if fin.any():
        g = _rc.reconstruct_3d(out1[fin], in1[fin], in2[fin], out2[fin],
                               sigma_hit, rng)
        sep[fin] = g['sep']; sep_outer[fin] = g['sep_outer']
        open_angle[fin] = g['open_angle']; dca[fin] = g['dca']
        collin[fin] = g['collin']; pointing[fin] = g['pointing']
        vtx_in[fin] = g['vtx_in']
        nf = int(fin.sum())
        true_hits = np.stack([out1[fin], in1[fin], in2[fin], out2[fin]], 1)
        smeared = np.stack([g['H_out1'], g['H_in1'], g['H_in2'], g['H_out2']], 1)
        tchi2, _, _, _ = _rc.timing_chi2_4hit(
            true_hits, decay_pos[fin], np.ones((nf, 4)), np.ones((nf, 4)),
            smeared, g['V_reco'], sigma_t, rng)
        timing[fin] = tchi2
    on_tracker = on_tracker & fin

    return dict(sep=sep, sep_outer=sep_outer, open_angle=open_angle, dca=dca,
                collin=collin, pointing=pointing, vtx_in=vtx_in,
                on_tracker=on_tracker, timing_chi2=timing,
                p_soft=np.asarray(p_soft, float))


def selection_mask(mc, p_cut=P_CUT, sep_min=SEP_MIN, sep_max=SEP_MAX,
                   dca_cut=DCA_CUT, theta_parallel=THETA_PARALLEL,
                   sep_out_max_parallel=SEP_OUT_MAX_PARALLEL,
                   collin_frac=COLLIN_FRAC, sep_out_gate=SEP_OUT_GATE,
                   point_tight_sep_in=POINT_TIGHT_SEP_IN,
                   sep_in_point_gate=SEP_IN_POINT_GATE,
                   point_global=POINT_GLOBAL,
                   apply_timing=True, chi2_timing_max=CHI2_TIMING_MAX):
    """Full PR #13 signal selection (mirrors higgs selection_mask on main)."""
    sep = mc['sep']; sep_outer = mc['sep_outer']; p_soft = mc['p_soft']
    dca = mc['dca']; open_angle = mc['open_angle']; collin = mc['collin']
    vtx_in = mc['vtx_in']; on_tracker = mc['on_tracker']; pointing = mc['pointing']

    with np.errstate(invalid='ignore'):
        m = (p_soft >= p_cut)
        m &= (sep >= sep_min) & (sep_outer >= sep_min)
        m &= (sep <= sep_max)
        m &= (dca <= dca_cut)
        is_parallel = open_angle < theta_parallel
        m &= (~is_parallel) | (sep_outer < sep_out_max_parallel)
        gated = sep_outer > sep_out_gate
        m &= (~gated) | (collin > collin_frac * DETECTOR_THICKNESS)
        gated_pt = sep < sep_in_point_gate
        m &= (~gated_pt) | (pointing < point_tight_sep_in)
        m &= pointing < point_global
        if apply_timing and 'timing_chi2' in mc:
            m &= mc['timing_chi2'] < chi2_timing_max
    m = m & vtx_in & on_tracker
    return m


# =========================================================================
# Per-event Monte Carlo and the U^2 scan
# =========================================================================
# The decay geometry (vertex, daughter directions, wall hits, reconstruction,
# selection) is independent of the mixing U^2, so the MC is built ONCE per mass
# point and reweighted to every U^2 by the decay-density factor, exactly as the
# higgs MC reweights a single uniform-sampled pass to each lifetime.

CMS_ORIGIN = np.zeros(3)


def build_event_mc(p4, direction, entry_d, exit_d, templates, n_samples, rng,
                   sigma_hit=HIT_RESOLUTION, sigma_t=SIGMA_T_DEFAULT,
                   origin=CMS_ORIGIN):
    """Sample decay vertices and reconstruct, for a batch of HNL four-vectors.

    For each of the ``n_events`` four-vectors, sample ``n_samples`` decay
    distances uniformly in ``[entry_d, exit_d]``, draw a FairShip decay template
    at each, boost it to the lab, take the best two charged tracks, and run the
    full reconstruction + selection. Decays with fewer than two charged tracks
    (the invisible fraction) fail, so the visible branching fraction is folded
    in here -- no separate BR_vis table is needed.

    Returns ``(d, passed)``, each ``(n_events, n_samples)``: the sampled decay
    distance and whether that decay passes the signal selection. Both are
    U^2-independent; ``scan_u2`` applies the lifetime weight.
    """
    origin = np.asarray(origin, float)
    p4 = np.asarray(p4, float)
    direction = np.asarray(direction, float)
    n_ev = len(entry_d)
    counts = templates['daughter_counts']
    off = np.concatenate([[0], np.cumsum(counts)])
    n_tmpl = len(counts)
    # Materialize the daughter arrays once. A compressed NpzFile re-decompresses
    # the whole array on every __getitem__, so slicing them inside the per-decay
    # loop below would otherwise decompress each array M times over.
    t_px = np.asarray(templates['px']); t_py = np.asarray(templates['py'])
    t_pz = np.asarray(templates['pz']); t_E = np.asarray(templates['energy'])
    t_charge = np.asarray(templates['charge'])
    t_stable = np.asarray(templates['stable'])

    d = rng.uniform(np.asarray(entry_d)[:, None], np.asarray(exit_d)[:, None],
                    size=(n_ev, n_samples))
    M = n_ev * n_samples
    vtx = origin[None, :] + d.reshape(M, 1) * np.repeat(direction, n_samples, axis=0)

    dir1 = np.empty((M, 3)); dir2 = np.empty((M, 3)); psoft = np.zeros(M)
    valid = np.zeros(M, dtype=bool)
    tmpl_idx = rng.integers(0, n_tmpl, size=M)
    for k in range(M):
        ti = tmpl_idx[k]
        s, e = off[ti], off[ti + 1]
        ev = k // n_samples
        lx, ly, lz, _ = boost_rest_to_lab(
            p4[ev], t_px[s:e], t_py[s:e], t_pz[s:e], t_E[s:e])
        bt = best_two_directions(lx, ly, lz, t_charge[s:e], t_stable[s:e])
        if bt is not None:
            dir1[k], dir2[k], psoft[k] = bt
            valid[k] = True

    passed = np.zeros(M, dtype=bool)
    if valid.any():
        mc = reconstruct_decays(vtx[valid], dir1[valid], dir2[valid],
                                psoft[valid], sigma_hit, sigma_t, rng)
        passed[valid] = selection_mask(mc)
    return d, passed.reshape(n_ev, n_samples)


def scan_u2(d, passed, path_len, weight, beta_gamma, ctau_u2_1,
            L_int_pb, u2_grid):
    """N_signal(U^2) by reweighting the once-built MC.

    For each event the decay-and-pass probability is the MC estimate of
    ``int (1/lam) e^{-x/lam} [pass] dx`` with ``lam = beta*gamma * ctau_u2_1/U^2``;
    ``N = L_int * U^2 * sum_events weight * P``. ``d``/``passed`` are
    ``(n_events, n_samples)`` from :func:`build_event_mc`.
    """
    n_ev, n_samples = d.shape
    pw = passed.astype(float)
    weight = np.asarray(weight, float)
    bg = np.asarray(beta_gamma, float)
    per_sample = (np.asarray(path_len, float) / n_samples)[:, None]
    N_grid = np.zeros(len(u2_grid))
    for iu, u2 in enumerate(u2_grid):
        lam = bg * ctau_u2_1 / u2                       # (n_ev,)
        inv = (1.0 / lam)[:, None]
        density = inv * np.exp(-d * inv)                # (n_ev, n_samples)
        P_ev = (per_sample * density * pw).sum(axis=1)  # (n_ev,)
        N_grid[iu] = L_int_pb * u2 * float(weight @ P_ev)
    return u2_grid, N_grid
