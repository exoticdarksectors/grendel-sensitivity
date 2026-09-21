"""Vectorized 2- and 3-body decay samplers for HNL production."""

import numpy as np


def _boost_to_lab(p4_rest, parent_E, parent_px, parent_py, parent_pz):
    """Boost rest-frame 4-vectors to the lab frame."""
    E_p = parent_E
    beta_x = parent_px / E_p
    beta_y = parent_py / E_p
    beta_z = parent_pz / E_p
    beta2 = beta_x**2 + beta_y**2 + beta_z**2
    beta2 = np.clip(beta2, 0, 1 - 1e-12)
    gamma = 1.0 / np.sqrt(1.0 - beta2)

    e_r = p4_rest[:, 0]
    px_r = p4_rest[:, 1]
    py_r = p4_rest[:, 2]
    pz_r = p4_rest[:, 3]

    bp = beta_x * px_r + beta_y * py_r + beta_z * pz_r

    at_rest = beta2 < 1e-30

    e_lab = gamma * (e_r + bp)
    px_lab = px_r.copy()
    py_lab = py_r.copy()
    pz_lab = pz_r.copy()

    moving = ~at_rest
    if np.any(moving):
        k = (gamma[moving] - 1.0) * bp[moving] / beta2[moving]
        px_lab[moving] += beta_x[moving] * k + gamma[moving] * beta_x[moving] * e_r[moving]
        py_lab[moving] += beta_y[moving] * k + gamma[moving] * beta_y[moving] * e_r[moving]
        pz_lab[moving] += beta_z[moving] * k + gamma[moving] * beta_z[moving] * e_r[moving]

    return np.column_stack([e_lab, px_lab, py_lab, pz_lab])


def decay_2body(parent_E, parent_px, parent_py, parent_pz, m_parent, m1, m2, rng=None):
    """Isotropic parent -> d1(m1) + d2(m2) decay."""
    if rng is None:
        rng = np.random.default_rng()

    N = len(parent_E)

    M = m_parent
    p_star = np.sqrt(
        np.maximum((M**2 - (m1 + m2)**2) * (M**2 - (m1 - m2)**2), 0.0)
    ) / (2.0 * M)

    cos_theta = rng.uniform(-1.0, 1.0, N)
    sin_theta = np.sqrt(1.0 - cos_theta**2)
    phi = rng.uniform(0.0, 2 * np.pi, N)

    px_star = p_star * sin_theta * np.cos(phi)
    py_star = p_star * sin_theta * np.sin(phi)
    pz_star = p_star * cos_theta

    E1_star_val = np.sqrt(p_star**2 + m1**2)
    E2_star_val = np.sqrt(p_star**2 + m2**2)
    E1_star = np.full(N, E1_star_val)
    E2_star = np.full(N, E2_star_val)

    d1_rest = np.column_stack([E1_star, px_star, py_star, pz_star])
    d2_rest = np.column_stack([E2_star, -px_star, -py_star, -pz_star])

    d1_lab = _boost_to_lab(d1_rest, parent_E, parent_px, parent_py, parent_pz)
    d2_lab = _boost_to_lab(d2_rest, parent_E, parent_px, parent_py, parent_pz)

    return d1_lab, d2_lab


def _sample_polar_cos(asym, n, rng):
    """Sample cos theta from pdf proportional to 1 + asym*cos(theta)."""
    u = rng.random(n)
    asym = float(np.clip(asym, -1.0, 1.0))
    if abs(asym) < 1e-12:
        return 2.0 * u - 1.0
    a = asym / 4.0
    b = 0.5
    c = 0.5 - asym / 4.0 - u
    disc = np.sqrt(np.maximum(b * b - 4.0 * a * c, 0.0))
    cos_theta = (-b + disc) / (2.0 * a)
    return np.clip(cos_theta, -1.0, 1.0)


def decay_2body_polarized(parent_E, parent_px, parent_py, parent_pz,
                          m_parent, m1, m2, asymmetry=0.0, rng=None):
    """Longitudinally polarized 2-body decay; d2 follows 1 + asym*cos(theta)."""
    if rng is None:
        rng = np.random.default_rng()

    N = len(parent_E)
    M = m_parent
    p_star = np.sqrt(
        np.maximum((M**2 - (m1 + m2)**2) * (M**2 - (m1 - m2)**2), 0.0)
    ) / (2.0 * M)

    p_mag = np.sqrt(parent_px**2 + parent_py**2 + parent_pz**2)
    safe = np.where(p_mag > 0, p_mag, 1.0)
    nx, ny, nz = parent_px / safe, parent_py / safe, parent_pz / safe
    at_rest = p_mag == 0
    nz = np.where(at_rest, 1.0, nz)
    nx = np.where(at_rest, 0.0, nx)
    ny = np.where(at_rest, 0.0, ny)

    cos_theta = _sample_polar_cos(asymmetry, N, rng)
    sin_theta = np.sqrt(np.maximum(1.0 - cos_theta**2, 0.0))
    phi = rng.uniform(0.0, 2 * np.pi, N)

    ref_z = np.abs(nz) < 0.9
    rx = np.where(ref_z, 0.0, 1.0)
    ry = np.zeros(N)
    rz = np.where(ref_z, 1.0, 0.0)
    e1x = ny * rz - nz * ry
    e1y = nz * rx - nx * rz
    e1z = nx * ry - ny * rx
    e1n = np.sqrt(e1x**2 + e1y**2 + e1z**2)
    e1n = np.where(e1n > 0, e1n, 1.0)
    e1x, e1y, e1z = e1x / e1n, e1y / e1n, e1z / e1n
    e2x = ny * e1z - nz * e1y
    e2y = nz * e1x - nx * e1z
    e2z = nx * e1y - ny * e1x

    dirx = cos_theta * nx + sin_theta * (np.cos(phi) * e1x + np.sin(phi) * e2x)
    diry = cos_theta * ny + sin_theta * (np.cos(phi) * e1y + np.sin(phi) * e2y)
    dirz = cos_theta * nz + sin_theta * (np.cos(phi) * e1z + np.sin(phi) * e2z)

    px2 = p_star * dirx
    py2 = p_star * diry
    pz2 = p_star * dirz
    E2_star = np.sqrt(p_star**2 + m2**2)
    E1_star = np.sqrt(p_star**2 + m1**2)

    d2_rest = np.column_stack([np.full(N, 0.0) + E2_star, px2, py2, pz2])
    d1_rest = np.column_stack([np.full(N, 0.0) + E1_star, -px2, -py2, -pz2])

    d1_lab = _boost_to_lab(d1_rest, parent_E, parent_px, parent_py, parent_pz)
    d2_lab = _boost_to_lab(d2_rest, parent_E, parent_px, parent_py, parent_pz)
    return d1_lab, d2_lab




def _eval_dbr(dbr_expr, locals_dict):
    """Vectorized eval of an HNLCalc differential-BR string."""
    if "mass" in locals_dict and "m3" not in locals_dict:
        locals_dict = dict(locals_dict, m3=locals_dict["mass"])
    return eval(dbr_expr, {"np": np, "__builtins__": {}}, locals_dict)


def _orthonormal_basis(nx, ny, nz):
    """For each unit vector, pick two perpendicular basis vectors."""
    n = len(nx)
    ref_z = np.abs(nz) < 0.9
    rx = np.where(ref_z, 0.0, 1.0)
    ry = np.zeros(n)
    rz = np.where(ref_z, 1.0, 0.0)
    e1x = ny * rz - nz * ry
    e1y = nz * rx - nx * rz
    e1z = nx * ry - ny * rx
    e1n = np.sqrt(e1x ** 2 + e1y ** 2 + e1z ** 2)
    e1n = np.where(e1n > 0, e1n, 1.0)
    e1x, e1y, e1z = e1x / e1n, e1y / e1n, e1z / e1n
    e2x = ny * e1z - nz * e1y
    e2y = nz * e1x - nx * e1z
    e2z = nx * e1y - ny * e1x
    return (e1x, e1y, e1z), (e2x, e2y, e2z)


def decay_3body_weighted_dq2dE(parent_E, parent_px, parent_py, parent_pz,
                               m_parent, m1, m2, m3, dbr_expr,
                               coupling=1.0, rng=None,
                               ceiling_grid=80, safety=1.2,
                               max_iter=10000):
    """Sample parent -> d1+d2+d3 with HNLCalc dBR/dq2/dE weighting."""
    if rng is None:
        rng = np.random.default_rng()

    N = len(parent_E)
    M = m_parent

    q2min = (m2 + m3) ** 2
    q2max = (M - m1) ** 2
    if q2max <= q2min:
        raise RuntimeError("No 3-body phase space (at or below threshold).")

    q2_grid = np.linspace(q2min + 1e-9, q2max - 1e-9, ceiling_grid)
    u_grid = np.linspace(1e-6, 1.0 - 1e-6, ceiling_grid)
    Q2, U = np.meshgrid(q2_grid, u_grid, indexing="ij")
    Q = np.sqrt(Q2)
    E2st = (Q2 - m2 ** 2 + m3 ** 2) / (2.0 * Q)
    E3st = (M ** 2 - Q2 - m1 ** 2) / (2.0 * Q)
    A = E2st + E3st
    B1 = np.sqrt(np.maximum(E2st ** 2 - m3 ** 2, 0.0))
    B2 = np.sqrt(np.maximum(E3st ** 2 - m1 ** 2, 0.0))
    m232min = A ** 2 - (B1 + B2) ** 2
    m232max = A ** 2 - (B1 - B2) ** 2
    ENmin_grid = (m232min + Q2 - m2 ** 2 - m1 ** 2) / (2.0 * M)
    ENmax_grid = (m232max + Q2 - m2 ** 2 - m1 ** 2) / (2.0 * M)
    E_grid = ENmin_grid + U * (ENmax_grid - ENmin_grid)
    w_grid = _eval_dbr(
        dbr_expr,
        {"q": Q, "q2": Q2, "energy": E_grid, "mass": m3, "coupling": coupling},
    )
    w_grid = np.where(np.isfinite(w_grid) & (w_grid > 0), w_grid, 0.0)
    span_grid = np.maximum(ENmax_grid - ENmin_grid, 0.0)
    ceiling = float((w_grid * span_grid).max()) * safety
    if not np.isfinite(ceiling) or ceiling <= 0:
        raise RuntimeError(
            f"dbr ceiling not positive (max={ceiling}). dbr_expr may be "
            f"degenerate at this (M, m1, m2, m3) = "
            f"({M}, {m1}, {m2}, {m3})."
        )

    q2_acc = np.empty(N)
    E_acc = np.empty(N)
    max_restarts = 5
    for _restart in range(max_restarts):
        accepted = np.zeros(N, dtype=bool)
        w_seen_max = 0.0
        for _ in range(max_iter):
            if accepted.all():
                break
            n_need = int((~accepted).sum())
            q2c = rng.uniform(q2min, q2max, n_need)
            qc = np.sqrt(q2c)
            E2 = (q2c - m2 ** 2 + m3 ** 2) / (2.0 * qc)
            E3 = (M ** 2 - q2c - m1 ** 2) / (2.0 * qc)
            A = E2 + E3
            B1 = np.sqrt(np.maximum(E2 ** 2 - m3 ** 2, 0.0))
            B2 = np.sqrt(np.maximum(E3 ** 2 - m1 ** 2, 0.0))
            mmin = A ** 2 - (B1 + B2) ** 2
            mmax = A ** 2 - (B1 - B2) ** 2
            ENmin = (mmin + q2c - m2 ** 2 - m1 ** 2) / (2.0 * M)
            ENmax = (mmax + q2c - m2 ** 2 - m1 ** 2) / (2.0 * M)
            uc = rng.random(n_need)
            Ec = ENmin + uc * (ENmax - ENmin)
            span = ENmax - ENmin
            wc = _eval_dbr(
                dbr_expr,
                {"q": qc, "q2": q2c, "energy": Ec, "mass": m3, "coupling": coupling},
            )
            wc = np.where(np.isfinite(wc) & (wc > 0), wc, 0.0)
            wspan = wc * span
            w_seen_max = max(w_seen_max, float(wspan.max()))
            u_acc = rng.uniform(0.0, ceiling, n_need)
            ok = u_acc < wspan
            idx = np.where(~accepted)[0]
            q2_acc[idx[ok]] = q2c[ok]
            E_acc[idx[ok]] = Ec[ok]
            accepted[idx[ok]] = True
        if w_seen_max > ceiling:
            # The grid scan under-estimated the true maximum, so the sample
            # is biased against the over-ceiling region: discard and restart
            # with the observed maximum.
            ceiling = w_seen_max * safety
            continue
        if not accepted.all():
            raise RuntimeError(
                f"dq2dE accept-reject failed to converge ({(~accepted).sum()} "
                f"of {N} unaccepted after {max_iter} iterations; ceiling may be "
                f"too loose)."
            )
        break
    else:
        raise RuntimeError(
            f"dq2dE ceiling still exceeded after {max_restarts} restarts; "
            f"dbr may be unbounded at this (M, m1, m2, m3) = "
            f"({M}, {m1}, {m2}, {m3})."
        )

    q_acc = np.sqrt(q2_acc)
    lam_d1 = (M ** 2 - (m1 + q_acc) ** 2) * (M ** 2 - (m1 - q_acc) ** 2)
    p_d1 = np.sqrt(np.maximum(lam_d1, 0.0)) / (2.0 * M)
    E_d1 = np.sqrt(p_d1 ** 2 + m1 ** 2)
    cos_t1 = rng.uniform(-1.0, 1.0, N)
    sin_t1 = np.sqrt(np.maximum(1.0 - cos_t1 ** 2, 0.0))
    phi1 = rng.uniform(0.0, 2 * np.pi, N)
    p_d1x = p_d1 * sin_t1 * np.cos(phi1)
    p_d1y = p_d1 * sin_t1 * np.sin(phi1)
    p_d1z = p_d1 * cos_t1

    E_sys23 = M - E_d1
    p_sys23x = -p_d1x
    p_sys23y = -p_d1y
    p_sys23z = -p_d1z
    gamma = E_sys23 / q_acc
    p_sys23_mag = p_d1  # magnitude, opposite direction
    beta = np.where(E_sys23 > 0, p_sys23_mag / E_sys23, 0.0)

    E_d3_star = (q2_acc + m3 ** 2 - m2 ** 2) / (2.0 * q_acc)
    p_d3_star = np.sqrt(np.maximum(E_d3_star ** 2 - m3 ** 2, 0.0))

    bg_p3 = beta * gamma * p_d3_star
    safe = bg_p3 > 1e-15
    cos_ts = np.where(safe, (E_acc - gamma * E_d3_star) / np.where(safe, bg_p3, 1.0), 0.0)
    cos_ts = np.clip(cos_ts, -1.0, 1.0)
    sin_ts = np.sqrt(np.maximum(1.0 - cos_ts ** 2, 0.0))
    phi_s = rng.uniform(0.0, 2 * np.pi, N)

    safe_d1 = p_d1 > 0
    nx = np.where(safe_d1, p_sys23x / np.where(safe_d1, p_sys23_mag, 1.0), 0.0)
    ny = np.where(safe_d1, p_sys23y / np.where(safe_d1, p_sys23_mag, 1.0), 0.0)
    nz = np.where(safe_d1, p_sys23z / np.where(safe_d1, p_sys23_mag, 1.0), 1.0)
    (e1x, e1y, e1z), (e2x, e2y, e2z) = _orthonormal_basis(nx, ny, nz)
    cphi = np.cos(phi_s)
    sphi = np.sin(phi_s)
    dx = cos_ts * nx + sin_ts * (cphi * e1x + sphi * e2x)
    dy = cos_ts * ny + sin_ts * (cphi * e1y + sphi * e2y)
    dz = cos_ts * nz + sin_ts * (cphi * e1z + sphi * e2z)

    p_d3x_star = p_d3_star * dx
    p_d3y_star = p_d3_star * dy
    p_d3z_star = p_d3_star * dz
    d3_in_sys23 = np.column_stack([
        np.full(N, 0.0) + E_d3_star,
        p_d3x_star, p_d3y_star, p_d3z_star,
    ])

    E_d2_star = np.sqrt(p_d3_star ** 2 + m2 ** 2)
    d2_in_sys23 = np.column_stack([
        np.full(N, 0.0) + E_d2_star,
        -p_d3x_star, -p_d3y_star, -p_d3z_star,
    ])

    d3_parent = _boost_to_lab(d3_in_sys23, E_sys23, p_sys23x, p_sys23y, p_sys23z)
    d2_parent = _boost_to_lab(d2_in_sys23, E_sys23, p_sys23x, p_sys23y, p_sys23z)
    d1_parent = np.column_stack([E_d1, p_d1x, p_d1y, p_d1z])

    d1_lab = _boost_to_lab(d1_parent, parent_E, parent_px, parent_py, parent_pz)
    d2_lab = _boost_to_lab(d2_parent, parent_E, parent_px, parent_py, parent_pz)
    d3_lab = _boost_to_lab(d3_parent, parent_E, parent_px, parent_py, parent_pz)

    return d1_lab, d2_lab, d3_lab


def decay_3body_weighted_dE(parent_E, parent_px, parent_py, parent_pz,
                            m_parent, m1, m2, m3, dbr_expr,
                            coupling=1.0, rng=None,
                            ceiling_grid=400, safety=1.2,
                            max_iter=10000):
    """Sample tau-leptonic 3-body decay with HNLCalc dBR/dE weighting."""
    if rng is None:
        rng = np.random.default_rng()

    N = len(parent_E)
    M = m_parent

    Emin = m3
    Emax = (M ** 2 + m3 ** 2 - (m1 + m2) ** 2) / (2.0 * M)
    if Emax <= Emin:
        raise RuntimeError("No 3-body phase space for dE sampler.")

    E_grid = np.linspace(Emin + 1e-9, Emax - 1e-9, ceiling_grid)
    w_grid = _eval_dbr(
        dbr_expr,
        {"energy": E_grid, "mass": m3, "coupling": coupling},
    )
    w_grid = np.where(np.isfinite(w_grid) & (w_grid > 0), w_grid, 0.0)
    ceiling = float(w_grid.max()) * safety
    if not np.isfinite(ceiling) or ceiling <= 0:
        raise RuntimeError(
            f"dE ceiling not positive (max={ceiling}); "
            f"(M, m1, m2, m3) = ({M}, {m1}, {m2}, {m3})."
        )

    E_acc = np.empty(N)
    max_restarts = 5
    for _restart in range(max_restarts):
        accepted = np.zeros(N, dtype=bool)
        w_seen_max = 0.0
        for _ in range(max_iter):
            if accepted.all():
                break
            n_need = int((~accepted).sum())
            Ec = rng.uniform(Emin, Emax, n_need)
            wc = _eval_dbr(
                dbr_expr,
                {"energy": Ec, "mass": m3, "coupling": coupling},
            )
            wc = np.where(np.isfinite(wc) & (wc > 0), wc, 0.0)
            w_seen_max = max(w_seen_max, float(wc.max()))
            u = rng.uniform(0.0, ceiling, n_need)
            ok = u < wc
            idx = np.where(~accepted)[0]
            E_acc[idx[ok]] = Ec[ok]
            accepted[idx[ok]] = True
        if w_seen_max > ceiling:
            # The grid scan under-estimated the true maximum, so the sample
            # is biased against the over-ceiling region: discard and restart
            # with the observed maximum.
            ceiling = w_seen_max * safety
            continue
        if not accepted.all():
            raise RuntimeError("dE accept-reject failed to converge.")
        break
    else:
        raise RuntimeError(
            f"dE ceiling still exceeded after {max_restarts} restarts; "
            f"(M, m1, m2, m3) = ({M}, {m1}, {m2}, {m3})."
        )

    p_d3 = np.sqrt(np.maximum(E_acc ** 2 - m3 ** 2, 0.0))
    cos_t3 = rng.uniform(-1.0, 1.0, N)
    sin_t3 = np.sqrt(np.maximum(1.0 - cos_t3 ** 2, 0.0))
    phi3 = rng.uniform(0.0, 2 * np.pi, N)
    p_d3x = p_d3 * sin_t3 * np.cos(phi3)
    p_d3y = p_d3 * sin_t3 * np.sin(phi3)
    p_d3z = p_d3 * cos_t3
    d3_parent = np.column_stack([E_acc, p_d3x, p_d3y, p_d3z])

    m12sq = M ** 2 + m3 ** 2 - 2.0 * M * E_acc
    m12 = np.sqrt(np.maximum(m12sq, 0.0))
    E_sys12 = M - E_acc
    lam12 = (m12sq - (m1 + m2) ** 2) * (m12sq - (m1 - m2) ** 2)
    p_star = np.sqrt(np.maximum(lam12, 0.0)) / (2.0 * np.maximum(m12, 1e-15))
    cos_a = rng.uniform(-1.0, 1.0, N)
    sin_a = np.sqrt(np.maximum(1.0 - cos_a ** 2, 0.0))
    phi_a = rng.uniform(0.0, 2 * np.pi, N)
    p1x = p_star * sin_a * np.cos(phi_a)
    p1y = p_star * sin_a * np.sin(phi_a)
    p1z = p_star * cos_a
    E1_star = np.sqrt(p_star ** 2 + m1 ** 2)
    E2_star = np.sqrt(p_star ** 2 + m2 ** 2)
    d1_in_sys12 = np.column_stack([E1_star, p1x, p1y, p1z])
    d2_in_sys12 = np.column_stack([E2_star, -p1x, -p1y, -p1z])
    d1_parent = _boost_to_lab(d1_in_sys12, E_sys12, -p_d3x, -p_d3y, -p_d3z)
    d2_parent = _boost_to_lab(d2_in_sys12, E_sys12, -p_d3x, -p_d3y, -p_d3z)

    d1_lab = _boost_to_lab(d1_parent, parent_E, parent_px, parent_py, parent_pz)
    d2_lab = _boost_to_lab(d2_parent, parent_E, parent_px, parent_py, parent_pz)
    d3_lab = _boost_to_lab(d3_parent, parent_E, parent_px, parent_py, parent_pz)

    return d1_lab, d2_lab, d3_lab


def _dalitz_q2_range(M, m1, m2, m3, s12):
    """(p2+p3)^2 range at fixed s12=(p1+p2)^2, mirroring HNLCalc dq2dm122."""
    sq = np.sqrt(s12)
    e2 = (s12 - m1 ** 2 + m2 ** 2) / (2.0 * sq)   # d2 energy in the (12) rest frame
    e3 = (M ** 2 - s12 - m3 ** 2) / (2.0 * sq)     # d3 energy in the (12) rest frame
    a = e2 + e3
    b1 = np.sqrt(np.maximum(e2 ** 2 - m2 ** 2, 0.0))
    b2 = np.sqrt(np.maximum(e3 ** 2 - m3 ** 2, 0.0))
    q2min = a ** 2 - (b1 + b2) ** 2
    q2max = a ** 2 - (b1 - b2) ** 2
    return q2min, q2max


def decay_3body_weighted_dq2dm122(parent_E, parent_px, parent_py, parent_pz,
                                  m_parent, m1, m2, m3, dbr_expr,
                                  coupling=1.0, rng=None,
                                  ceiling_grid=80, safety=1.2, max_iter=10000):
    """Sample parent -> d1+d2+d3 with HNLCalc dBR/dq2/dm12^2 weighting.

    This mirrors ``HNLCalc.integrate_3body_br_3body_dq2dm122`` (used by the
    baryon channels): the differential rate is parametrized in
    ``m12sq = (p1+p2)^2`` and ``q2 = (p2+p3)^2`` rather than in the HNL energy.
    ``d3`` is the HNL (mass ``m3``). The dbr string references ``q2``,
    ``m12sq``, ``mass`` and ``coupling``. Returns ``(d1_lab, d2_lab, d3_lab)``.
    """
    if rng is None:
        rng = np.random.default_rng()

    N = len(parent_E)
    M = m_parent

    s12min = (m1 + m2) ** 2
    s12max = (M - m3) ** 2
    if s12max <= s12min:
        raise RuntimeError("No 3-body phase space (dq2dm122, at or below threshold).")

    # Accept-reject ceiling on a grid over (s12, q2-fraction).
    s12_grid = np.linspace(s12min + 1e-9, s12max - 1e-9, ceiling_grid)
    u_grid = np.linspace(1e-6, 1.0 - 1e-6, ceiling_grid)
    S12, U = np.meshgrid(s12_grid, u_grid, indexing="ij")
    q2min_g, q2max_g = _dalitz_q2_range(M, m1, m2, m3, S12)
    Q2g = q2min_g + U * (q2max_g - q2min_g)
    span_g = np.maximum(q2max_g - q2min_g, 0.0)
    wg = _eval_dbr(dbr_expr, {"q2": Q2g, "m12sq": S12, "mass": m3, "coupling": coupling})
    wg = np.where(np.isfinite(wg) & (wg > 0), wg, 0.0)
    ceiling = float((wg * span_g).max()) * safety
    if not np.isfinite(ceiling) or ceiling <= 0:
        raise RuntimeError(
            f"dq2dm122 ceiling not positive (max={ceiling}); "
            f"(M, m1, m2, m3) = ({M}, {m1}, {m2}, {m3})."
        )

    s12_acc = np.empty(N)
    q2_acc = np.empty(N)
    max_restarts = 5
    for _restart in range(max_restarts):
        accepted = np.zeros(N, dtype=bool)
        w_seen_max = 0.0
        for _ in range(max_iter):
            if accepted.all():
                break
            n_need = int((~accepted).sum())
            s12c = rng.uniform(s12min, s12max, n_need)
            q2min_c, q2max_c = _dalitz_q2_range(M, m1, m2, m3, s12c)
            uc = rng.random(n_need)
            q2c = q2min_c + uc * (q2max_c - q2min_c)
            span = q2max_c - q2min_c
            wc = _eval_dbr(dbr_expr, {"q2": q2c, "m12sq": s12c, "mass": m3, "coupling": coupling})
            wc = np.where(np.isfinite(wc) & (wc > 0), wc, 0.0)
            wspan = wc * span
            w_seen_max = max(w_seen_max, float(wspan.max()))
            ua = rng.uniform(0.0, ceiling, n_need)
            ok = ua < wspan
            idx = np.where(~accepted)[0]
            s12_acc[idx[ok]] = s12c[ok]
            q2_acc[idx[ok]] = q2c[ok]
            accepted[idx[ok]] = True
        if w_seen_max > ceiling:
            # The grid scan under-estimated the true maximum, so the sample
            # is biased against the over-ceiling region: discard and restart
            # with the observed maximum.
            ceiling = w_seen_max * safety
            continue
        if not accepted.all():
            raise RuntimeError(
                f"dq2dm122 accept-reject failed to converge "
                f"({(~accepted).sum()} of {N} unaccepted after {max_iter} iterations)."
            )
        break
    else:
        raise RuntimeError(
            f"dq2dm122 ceiling still exceeded after {max_restarts} restarts; "
            f"(M, m1, m2, m3) = ({M}, {m1}, {m2}, {m3})."
        )

    # Dalitz reconstruction in the parent rest frame from (s12, s23=q2).
    s12 = s12_acc
    s23 = q2_acc
    s13 = M ** 2 + m1 ** 2 + m2 ** 2 + m3 ** 2 - s12 - s23

    E1 = (M ** 2 + m1 ** 2 - s23) / (2.0 * M)
    E2 = (M ** 2 + m2 ** 2 - s13) / (2.0 * M)
    E3 = (M ** 2 + m3 ** 2 - s12) / (2.0 * M)
    p1 = np.sqrt(np.maximum(E1 ** 2 - m1 ** 2, 0.0))
    p3 = np.sqrt(np.maximum(E3 ** 2 - m3 ** 2, 0.0))

    # Angle between p1 and p3 from s13 = m1^2 + m3^2 + 2(E1 E3 - p1 p3 cos13).
    denom = p1 * p3
    safe = denom > 1e-15
    cos13 = np.where(
        safe,
        (E1 * E3 - (s13 - m1 ** 2 - m3 ** 2) / 2.0) / np.where(safe, denom, 1.0),
        0.0,
    )
    cos13 = np.clip(cos13, -1.0, 1.0)
    sin13 = np.sqrt(np.maximum(1.0 - cos13 ** 2, 0.0))

    # Isotropic orientation: random direction for p1, random azimuth for p3.
    cos_t = rng.uniform(-1.0, 1.0, N)
    sin_t = np.sqrt(np.maximum(1.0 - cos_t ** 2, 0.0))
    phi = rng.uniform(0.0, 2 * np.pi, N)
    n1x = sin_t * np.cos(phi)
    n1y = sin_t * np.sin(phi)
    n1z = cos_t
    (e1x, e1y, e1z), (e2x, e2y, e2z) = _orthonormal_basis(n1x, n1y, n1z)
    psi = rng.uniform(0.0, 2 * np.pi, N)
    cpsi, spsi = np.cos(psi), np.sin(psi)
    d3x = cos13 * n1x + sin13 * (cpsi * e1x + spsi * e2x)
    d3y = cos13 * n1y + sin13 * (cpsi * e1y + spsi * e2y)
    d3z = cos13 * n1z + sin13 * (cpsi * e1z + spsi * e2z)

    p1x, p1y, p1z = p1 * n1x, p1 * n1y, p1 * n1z
    p3x, p3y, p3z = p3 * d3x, p3 * d3y, p3 * d3z
    p2x, p2y, p2z = -(p1x + p3x), -(p1y + p3y), -(p1z + p3z)

    d1_rest = np.column_stack([E1, p1x, p1y, p1z])
    d2_rest = np.column_stack([E2, p2x, p2y, p2z])
    d3_rest = np.column_stack([E3, p3x, p3y, p3z])

    d1_lab = _boost_to_lab(d1_rest, parent_E, parent_px, parent_py, parent_pz)
    d2_lab = _boost_to_lab(d2_rest, parent_E, parent_px, parent_py, parent_pz)
    d3_lab = _boost_to_lab(d3_rest, parent_E, parent_px, parent_py, parent_pz)

    return d1_lab, d2_lab, d3_lab
