"""
Shared 4-hit reconstruction (single source of truth).

Both the signal (LLP->e+e-) and the cosmic-decay background reconstruct a
two-track vertex from four 3D hits: each track has an inner hit (closer to the
vertex) and an outer hit (farther). This module computes every geometric
observable (separations, opening angle, DCA, collinearity, PCA vertex, pointing)
identically from REAL 3D hit positions, so signal and background can never
disagree on a variable definition.

Convention: pass (A_out, A_in, B_in, B_out) where *_out is the hit farther from
the vertex on each track. Directions (out - in) then point away from the vertex,
so the pointing bisector is the sum of the two outgoing unit directions.
"""
import numpy as np
from .grendel_geometry import (classify_points_with_basis, DETECTOR_THICKNESS,
                              points_in_fiducial, points_on_tracker,
                              SPEED_OF_LIGHT)

IP = np.array([0.0, 0.0, 0.0])

# Margin used by wall_inner_outer to push the radial-inward projection safely
# PAST the fiducial face into the air gap (so the in-fiducial test is unambiguous).
# Derived from the geometry module's shell thickness rather than a hard-coded
# 10 cm, so it scales with the detector if DETECTOR_THICKNESS changes; it only
# needs to stay below the local air-gap depth (the gap is >> 0.4 L everywhere).
_FIDUCIAL_PROJ_MARGIN = 0.4 * DETECTOR_THICKNESS   # ~10 cm at L = 24 cm

# Timing-consistency cut (single source for signal + cosmic background).
SIGMA_T_DEFAULT = 0.5e-9      # per-hit timing resolution (s)
CHI2_TIMING_MAX = 6.0         # chi2 (ndof = 3) cut; chi2(3) CDF(6) ~ 0.89 = signal timing eff


def timing_chi2_4hit(true_hits, V_true, beta, sign, smeared_hits, V_reco,
                     sigma_t, rng):
    """Timing-consistency chi2 (ndof = 3) under the hypothesis that all four
    hits come from particles leaving the *reconstructed* vertex at c.

    Identical model for signal and the cosmic background:
      - true arrival time of each hit (t = 0 at the true vertex) is
        sign * r_true / (beta c), with sign = +1 for a hit AFTER the vertex
        (outgoing track) and -1 BEFORE it (incoming track, e.g. the cosmic
        muon); r_true = true-hit-to-true-vertex distance, beta = particle speed.
      - times are smeared by sigma_t.
      - predicted time t_i = t0 + R_i / c with R_i = smeared-hit-to-reco-vertex
        distance; the single free t0 (= mean residual) is fit out (ndof = 3).

    true_hits, smeared_hits: (n, 4, 3).  beta, sign: (n, 4).
    Returns (chi2, rms, t_meas, R) so callers can also build per-track tests.
    """
    c = SPEED_OF_LIGHT
    r_true = np.linalg.norm(true_hits - V_true[:, None, :], axis=2)
    t_true = sign * r_true / (beta * c)
    t_meas = t_true + rng.normal(0, sigma_t, t_true.shape)
    R = np.linalg.norm(smeared_hits - V_reco[:, None, :], axis=2)
    x = t_meas - R / c
    resid = x - x.mean(axis=1, keepdims=True)
    chi2 = (resid**2).sum(axis=1) / sigma_t**2
    rms = np.sqrt((resid**2).mean(axis=1))
    return chi2, rms, t_meas, R


def collinearity_4hit(pts):
    """RMS perpendicular residual of 4 hits (N,4,3) to their best-fit line."""
    pts = np.asarray(pts, float)
    ctr = pts.mean(axis=1, keepdims=True)
    _, S, _ = np.linalg.svd(pts - ctr, full_matrices=False)
    return np.sqrt((S[:, 1]**2 + S[:, 2]**2) / 4)


def _transverse_basis(d):
    helper = np.tile(np.array([0., 1., 0.]), (len(d), 1))
    par = np.abs((d * helper).sum(axis=1)) > 0.95
    helper[par] = np.array([1., 0., 0.])
    e1 = np.cross(d, helper)
    e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
    e2 = np.cross(d, e1)
    return e1, e2


def _line_closest_approach(a1, d1, a2, d2):
    w0 = a1 - a2
    a = np.einsum('ij,ij->i', d1, d1)
    b = np.einsum('ij,ij->i', d1, d2)
    c = np.einsum('ij,ij->i', d2, d2)
    dd = np.einsum('ij,ij->i', d1, w0)
    e = np.einsum('ij,ij->i', d2, w0)
    den = a * c - b * b
    den = np.where(np.abs(den) < 1e-30, 1e-30, den)
    t = (b * e - c * dd) / den
    s = (a * e - b * dd) / den
    return 0.5 * ((a1 + t[:, None] * d1) + (a2 + s[:, None] * d2))


def reconstruct_3d(A_out, A_in, B_in, B_out, sigma_hit, rng):
    """
    Geometric reconstruction from 4 real 3D hits. *_out = farther from vertex.
    Returns a dict with signal-convention keys: sep, sep_outer, open_angle, dca,
    collin, pointing, vtx_in, V_reco, and the smeared hits H_* + directions.
    """
    A_out, A_in = np.asarray(A_out, float), np.asarray(A_in, float)
    B_in, B_out = np.asarray(B_in, float), np.asarray(B_out, float)
    n = len(A_in)

    d1t = A_out - A_in
    d2t = B_out - B_in
    d1t /= np.linalg.norm(d1t, axis=1, keepdims=True)
    d2t /= np.linalg.norm(d2t, axis=1, keepdims=True)

    def smear(P, d):
        if sigma_hit <= 0:
            return P.copy()
        e1, e2 = _transverse_basis(d)
        return (P + rng.normal(0, sigma_hit, (n, 1)) * e1
                + rng.normal(0, sigma_hit, (n, 1)) * e2)

    Ho1, Hi1 = smear(A_out, d1t), smear(A_in, d1t)
    Hi2, Ho2 = smear(B_in, d2t), smear(B_out, d2t)

    d1 = Ho1 - Hi1
    d2 = Ho2 - Hi2
    m1 = np.linalg.norm(d1, axis=1)
    m2 = np.linalg.norm(d2, axis=1)

    sep = np.linalg.norm(Hi1 - Hi2, axis=1)
    sep_outer = np.linalg.norm(Ho1 - Ho2, axis=1)
    open_angle = np.arccos(np.clip((d1 * d2).sum(axis=1) / (m1 * m2), -1, 1))

    nvec = np.cross(d1, d2)
    nm = np.linalg.norm(nvec, axis=1)
    w = Hi1 - Hi2
    dca = np.where(nm > 1e-20,
                   np.abs((w * nvec).sum(axis=1)) / nm,
                   np.linalg.norm(np.cross(w, d1 / m1[:, None]), axis=1))

    collin = collinearity_4hit(np.stack([Ho1, Hi1, Hi2, Ho2], axis=1))

    V_reco = _line_closest_approach(Hi1, d1, Hi2, d2)
    vtx_in = points_in_fiducial(V_reco)

    u1 = d1 / m1[:, None]
    u2 = d2 / m2[:, None]
    bis = u1 + u2
    bm = np.linalg.norm(bis, axis=1)
    flight = V_reco - IP
    flight /= np.linalg.norm(flight, axis=1, keepdims=True)
    cosp = np.where(bm > 1e-12, (bis * flight).sum(axis=1) / bm, 0.0)
    pointing = np.arccos(np.clip(cosp, -1, 1))

    return dict(sep=sep, sep_outer=sep_outer, open_angle=open_angle, dca=dca,
                collin=collin, pointing=pointing, vtx_in=vtx_in, V_reco=V_reco,
                H_out1=Ho1, H_in1=Hi1, H_in2=Hi2, H_out2=Ho2, d1=d1, d2=d2)


def wall_inner_outer(exit_pt, direction, L=DETECTOR_THICKNESS):
    """
    Given a track's wall-crossing point (inner hit) and its 3D direction,
    return (inner, outer) real 3D hits, the outer being L farther along the
    track measured radially (the layer-2 hit).

    The along-track stub is L / |d.n_hat| (radial layer spacing L divided by the
    track's normal cosine), so an oblique crossing has a long stub. That stub is
    only physical while the track stays inside the instrumented shell of the
    cavern: a grazing track exits the cavern before reaching the outer layer and
    makes only ONE hit. We enforce this by projecting the candidate outer hit
    radially back onto the fiducial face and requiring that point inside the
    fiducial volume; tracks failing it (grazing / would-exit, and the
    ill-defined-normal case) return NaN for the outer hit -> unreconstructable.
    The outer hit must ALSO land on a TRACKER wall (points_on_tracker): an oblique
    track whose inner hit is on the tracker but whose layer-2 crossing walks over
    the tracker/scintillator boundary does not make two tracker hits -- counting
    it would be optimistic for signal (one real tracker hit) and (anti-)conservative
    for background (the real event fires the scintillator veto). This keeps the
    4-hit geometry physical for BOTH signal and the cosmic background (single
    source) and prevents grazing tracks from inflating sep_outer / collinearity.
    """
    exit_pt = np.asarray(exit_pt, float)
    direction = np.asarray(direction, float)
    theta, _, _, right, up = classify_points_with_basis(exit_pt)
    n_hat = np.cos(theta)[:, None] * right + np.sin(theta)[:, None] * up
    d_dot_n = np.einsum('ij,ij->i', direction, n_hat)
    d_dot_n = np.where(np.abs(d_dot_n) < 1e-6, np.nan, d_dot_n)
    outer = exit_pt + (L / d_dot_n)[:, None] * direction

    # Physical bound: the outer hit must sit on the cavern shell (project it
    # radially inward, just past the fiducial face, and require it in-fiducial)
    # AND land on a tracker wall (not the scintillator complement).
    finite = np.isfinite(outer[:, 0])
    inside = np.zeros(len(outer), dtype=bool)
    if finite.any():
        th_o, _, _, right_o, up_o = classify_points_with_basis(outer[finite])
        n_o = np.cos(th_o)[:, None] * right_o + np.sin(th_o)[:, None] * up_o
        # Go inward by L (the outer hit's radial offset from the fiducial face)
        # plus a geometry-derived margin to land safely inside the air gap.
        face_pt = outer[finite] - (L + _FIDUCIAL_PROJ_MARGIN) * n_o
        inside[finite] = (points_in_fiducial(face_pt)
                          & points_on_tracker(outer[finite]))
    outer = np.where(inside[:, None], outer, np.nan)
    return exit_pt, outer
