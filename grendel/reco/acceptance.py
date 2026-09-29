"""Signal acceptance from cached rest-frame decay templates."""
from __future__ import annotations

import os

import numpy as np

from ..constants import CMS_ORIGIN
from ..geometry.grendel_geometry import (
    mesh_fiducial, points_on_tracker, DETECTOR_THICKNESS,
)
from ..geometry import reco_common as _rc
from ..geometry.reco_common import SIGMA_T_DEFAULT, CHI2_TIMING_MAX

P_CUT = float(os.environ.get("GRENDEL_TRACK_P_CUT", "0.100"))
SEP_MIN = 0.01
SEP_MAX = 10.0
DCA_CUT = 0.1
THETA_PARALLEL = 0.100
SEP_OUT_MAX_PARALLEL = 0.30
COLLIN_FRAC = 0.40
SEP_OUT_GATE = 0.16
POINT_TIGHT_SEP_IN = 0.050
SEP_IN_POINT_GATE = 0.10
POINT_GLOBAL = 0.8
HIT_RESOLUTION = 0.003


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
    """General Lorentz boost of rest-frame daughter 3-momenta into the lab, along the parent
    3-momentum (replicates fairship_decay.boost_decay_to_lab on flat arrays)."""
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


def best_two_directions(px, py, pz, E, charge, stable, p_cut=P_CUT):
    """From one decay's lab daughters, pick the two highest-momentum charged stable tracks above
    p_cut."""
    p = np.sqrt(px**2 + py**2 + pz**2)
    sel = (np.abs(charge) > 0.5) & stable.astype(bool) & (p > p_cut)
    if int(sel.sum()) < 2:
        return None
    idx = np.where(sel)[0]
    order = idx[np.argsort(p[idx])[::-1]]
    i1, i2 = order[0], order[1]
    d1 = np.array([px[i1], py[i1], pz[i1]]) / p[i1]
    d2 = np.array([px[i2], py[i2], pz[i2]]) / p[i2]
    E = np.asarray(E, float)
    b1 = float(p[i1] / max(E[i1], p[i1]))
    b2 = float(p[i2] / max(E[i2], p[i2]))
    return d1, d2, float(min(p[i1], p[i2])), b1, b2


def reconstruct_decays(decay_pos, dir1, dir2, p_soft, sigma_hit, sigma_t, rng,
                       beta1=None, beta2=None, mesh=mesh_fiducial):
    """Build the selection ``mc`` dict for N decays from their vertex + the two daughter directions
    (the shared geometry->reco core)."""
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
        b1 = np.ones(n) if beta1 is None else np.asarray(beta1, float)
        b2 = np.ones(n) if beta2 is None else np.asarray(beta2, float)
        beta4 = np.column_stack([b1[fin], b1[fin], b2[fin], b2[fin]])
        tchi2, _, _, _ = _rc.timing_chi2_4hit(
            true_hits, decay_pos[fin], beta4, np.ones((nf, 4)),
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
    """Full signal selection: gate, pointing, collinearity and timing cuts."""
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


def scatter_mc(mc, valid, n_total):
    """Full-length selection arrays for the sequential cutflow."""
    full = {
        'sep': np.full(n_total, np.nan),
        'sep_outer': np.full(n_total, np.nan),
        'open_angle': np.full(n_total, np.nan),
        'dca': np.full(n_total, np.nan),
        'collin': np.full(n_total, np.nan),
        'pointing': np.full(n_total, np.nan),
        'vtx_in': np.zeros(n_total, dtype=bool),
        'on_tracker': np.zeros(n_total, dtype=bool),
        'timing_chi2': np.full(n_total, np.inf),
        'p_soft': np.zeros(n_total),
        'visible': np.asarray(valid, dtype=bool).copy(),
    }
    if mc is not None:
        for key in ('sep', 'sep_outer', 'open_angle', 'dca', 'collin',
                    'pointing', 'vtx_in', 'on_tracker', 'timing_chi2',
                    'p_soft'):
            full[key][valid] = mc[key]
    return full


def build_cutflow_mc(mc, weights, p_cut=P_CUT, sep_min=SEP_MIN,
                     sep_max=SEP_MAX, dca_cut=DCA_CUT,
                     theta_parallel=THETA_PARALLEL,
                     sep_out_max_parallel=SEP_OUT_MAX_PARALLEL,
                     collin_frac=COLLIN_FRAC, sep_out_gate=SEP_OUT_GATE,
                     point_tight_sep_in=POINT_TIGHT_SEP_IN,
                     sep_in_point_gate=SEP_IN_POINT_GATE,
                     point_global=POINT_GLOBAL,
                     chi2_timing_max=CHI2_TIMING_MAX):
    """Sequential signal cutflow from :func:`scatter_mc` full-length arrays."""
    weights = np.asarray(weights, float)
    total_w = weights.sum()
    if total_w <= 0:
        return []
    rows = []
    mask = np.ones(len(weights), dtype=bool)

    def add_row(name, new_mask, prev_mask):
        w = weights[new_mask].sum()
        prev_w = weights[prev_mask].sum()
        rows.append({'cut': name, 'weighted_yield': w,
                     'efficiency': w / total_w,
                     'marginal_efficiency': w / prev_w if prev_w > 0 else 0.0})

    add_row('Decay in fiducial volume', mask, mask)
    with np.errstate(invalid='ignore'):
        steps = [
            (f'>= 2 charged daughters (p > {p_cut*1000:.0f} MeV)',
             mc['visible']),
            ('Both daughters on tracker', mc['on_tracker']),
            (f'sep_in & sep_out > {sep_min*100:.0f} cm',
             (mc['sep'] >= sep_min) & (mc['sep_outer'] >= sep_min)),
            (f'sep_in < {sep_max:.0f} m', mc['sep'] <= sep_max),
            (f'DCA < {dca_cut*100:.0f} cm', mc['dca'] <= dca_cut),
            (f'sep_out < {sep_out_max_parallel*100:.0f} cm if theta < '
             f'{theta_parallel*1000:.0f} mrad',
             ~(mc['open_angle'] < theta_parallel)
             | (mc['sep_outer'] < sep_out_max_parallel)),
            (f'collin > {collin_frac*DETECTOR_THICKNESS*1000:.0f} mm if '
             f'sep_out > {sep_out_gate*100:.0f} cm',
             ~(mc['sep_outer'] > sep_out_gate)
             | (mc['collin'] > collin_frac * DETECTOR_THICKNESS)),
            ('Vertex in fiducial (PCA)', mc['vtx_in']),
            (f'pointing < {point_tight_sep_in*1000:.0f} mrad if sep_in < '
             f'{sep_in_point_gate*100:.0f} cm',
             ~(mc['sep'] < sep_in_point_gate)
             | (mc['pointing'] < point_tight_sep_in)),
            (f'pointing < {point_global*1000:.0f} mrad',
             mc['pointing'] < point_global),
            (f'timing chi2 < {chi2_timing_max:.0f}',
             mc['timing_chi2'] < chi2_timing_max),
        ]
        for name, ok in steps:
            prev = mask
            mask = mask & ok
            add_row(name, mask, prev)
    return rows


def build_event_mc(p4, direction, entry_d, exit_d, templates, n_samples, rng,
                   sigma_hit=HIT_RESOLUTION, sigma_t=SIGMA_T_DEFAULT,
                   origin=CMS_ORIGIN, return_mc=False):
    """Sample decay vertices and reconstruct, for a batch of HNL four-vectors."""
    origin = np.asarray(origin, float)
    p4 = np.asarray(p4, float)
    direction = np.asarray(direction, float)
    n_ev = len(entry_d)
    counts = templates['daughter_counts']
    off = np.concatenate([[0], np.cumsum(counts)])
    n_tmpl = len(counts)
    t_px = np.asarray(templates['px']); t_py = np.asarray(templates['py'])
    t_pz = np.asarray(templates['pz']); t_E = np.asarray(templates['energy'])
    t_charge = np.asarray(templates['charge'])
    t_stable = np.asarray(templates['stable'])

    d = rng.uniform(np.asarray(entry_d)[:, None], np.asarray(exit_d)[:, None],
                    size=(n_ev, n_samples))
    M = n_ev * n_samples
    vtx = origin[None, :] + d.reshape(M, 1) * np.repeat(direction, n_samples, axis=0)

    dir1 = np.empty((M, 3)); dir2 = np.empty((M, 3)); psoft = np.zeros(M)
    beta1 = np.ones(M); beta2 = np.ones(M)
    valid = np.zeros(M, dtype=bool)
    tmpl_idx = rng.integers(0, n_tmpl, size=M)
    for k in range(M):
        ti = tmpl_idx[k]
        s, e = off[ti], off[ti + 1]
        ev = k // n_samples
        lx, ly, lz, lE = boost_rest_to_lab(
            p4[ev], t_px[s:e], t_py[s:e], t_pz[s:e], t_E[s:e])
        bt = best_two_directions(lx, ly, lz, lE, t_charge[s:e], t_stable[s:e])
        if bt is not None:
            dir1[k], dir2[k], psoft[k], beta1[k], beta2[k] = bt
            valid[k] = True

    passed = np.zeros(M, dtype=bool)
    mc = None
    if valid.any():
        mc = reconstruct_decays(vtx[valid], dir1[valid], dir2[valid],
                                psoft[valid], sigma_hit, sigma_t, rng,
                                beta1=beta1[valid], beta2=beta2[valid])
        passed[valid] = selection_mask(mc)
    out = (d, passed.reshape(n_ev, n_samples), tmpl_idx.reshape(n_ev, n_samples))
    if return_mc:
        return (*out, scatter_mc(mc, valid, M))
    return out


def classify_template_modes(templates):
    """Per-template hadronic flag (length n_templates): True if the decay's daughters include any
    hadron (|pdg| > 100), else leptonic (leptons / neutrinos / photon only)."""
    counts = np.asarray(templates["daughter_counts"])
    pdg = np.abs(np.asarray(templates["pdg"]))
    off = np.concatenate([[0], np.cumsum(counts)])
    is_had = np.zeros(len(counts), dtype=bool)
    for i in range(len(counts)):
        seg = pdg[off[i]:off[i + 1]]
        is_had[i] = bool(seg.size and seg.max() > 100)
    return is_had


def scan_u2(d, passed, path_len, weight, beta_gamma, ctau_u2_1,
            L_int_pb, u2_grid, sample_w=None):
    """N_signal(U^2) by reweighting the once-built MC."""
    n_ev, n_samples = d.shape
    pw = passed.astype(float)
    if sample_w is not None:
        pw = pw * np.asarray(sample_w, float)
    weight = np.asarray(weight, float)
    bg = np.asarray(beta_gamma, float)
    per_sample = (np.asarray(path_len, float) / n_samples)[:, None]
    N_grid = np.zeros(len(u2_grid))
    for iu, u2 in enumerate(u2_grid):
        lam = bg * ctau_u2_1 / u2
        inv = (1.0 / lam)[:, None]
        density = inv * np.exp(-d * inv)
        P_ev = (per_sample * density * pw).sum(axis=1)
        N_grid[iu] = L_int_pb * u2 * float(weight @ P_ev)
    return u2_grid, N_grid


def signal_contribution_diagnostics(
    d,
    passed,
    path_len,
    weight,
    beta_gamma,
    ctau_u2_1,
    u2,
    sample_w=None,
):
    """Effective statistics of the weighted signal estimator at one coupling."""
    d = np.asarray(d, dtype=float)
    passed = np.asarray(passed, dtype=float)
    if d.ndim != 2 or passed.shape != d.shape:
        raise ValueError("d and passed must be equal-shape 2D arrays")
    n_events, n_samples = d.shape
    path_len = np.asarray(path_len, dtype=float)
    weight = np.asarray(weight, dtype=float)
    beta_gamma = np.asarray(beta_gamma, dtype=float)
    if any(len(array) != n_events for array in (path_len, weight, beta_gamma)):
        raise ValueError("event arrays must match the first dimension of d")
    if not np.isfinite(u2) or u2 <= 0.0:
        raise ValueError("u2 must be finite and positive")

    lifetime = beta_gamma * ctau_u2_1 / u2
    contribution = (
        weight[:, None]
        * (path_len / n_samples)[:, None]
        * np.exp(-d / lifetime[:, None])
        / lifetime[:, None]
        * passed
    )
    if sample_w is not None:
        contribution *= np.asarray(sample_w, dtype=float)
    total = float(contribution.sum())
    event_contribution = contribution.sum(axis=1)

    def ess(values):
        denominator = float(np.square(values).sum())
        return total**2 / denominator if denominator > 0.0 else 0.0

    return {
        "sample_ess": ess(contribution),
        "event_ess": ess(event_contribution),
        "max_event_fraction": (
            float(event_contribution.max()) / total if total > 0.0 else np.nan
        ),
        "nonzero_samples": int(np.count_nonzero(contribution)),
        "nonzero_events": int(np.count_nonzero(event_contribution)),
    }
