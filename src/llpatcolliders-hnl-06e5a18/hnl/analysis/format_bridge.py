"""
analysis/format_bridge.py

Convert HNL 4-vector CSVs (weight, E, px, py, pz) to geometry-ready arrays
(eta, phi, p, mass, beta, gamma, beta_gamma) for use with the GRENDEL
ray-casting code.
"""

import numpy as np
from pathlib import Path


def load_combined_csv(csv_path, m_N):
    """
    Load a combined HNL 4-vector CSV and compute derived kinematics.

    Parameters
    ----------
    csv_path : Path or str
        Headerless CSV with columns: weight, E, px, py, pz
    m_N : float
        HNL mass in GeV (from the mass grid).

    Returns
    -------
    dict with keys:
        weight     : (N,) production weight in pb at U²=1
        E, px, py, pz : (N,) 4-momentum components in GeV
        p          : (N,) total momentum |p|
        pt         : (N,) transverse momentum
        eta        : (N,) pseudorapidity
        phi        : (N,) azimuthal angle [−π, π)
        gamma      : (N,) Lorentz gamma factor
        beta       : (N,) velocity β = p/E
        beta_gamma : (N,) p / m_N
        mass       : (N,) filled with m_N
    """
    csv_path = Path(csv_path)
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return _empty_dict()

    data = np.loadtxt(csv_path, delimiter=",")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if len(data) == 0:
        return _empty_dict()

    weight = data[:, 0]
    E = data[:, 1]
    px = data[:, 2]
    py = data[:, 3]
    pz = data[:, 4]

    p = np.sqrt(px**2 + py**2 + pz**2)
    pt = np.sqrt(px**2 + py**2)

    # Pseudorapidity via polar angle
    theta = np.arctan2(pt, pz)
    theta = np.clip(theta, 1e-10, np.pi - 1e-10)
    eta = -np.log(np.tan(theta / 2.0))

    phi = np.arctan2(py, px)

    # Relativistic factors
    energy = np.maximum(E, m_N)  # guard against numerical noise
    gamma = energy / m_N
    beta = p / energy
    beta_gamma = p / m_N

    return {
        "weight": weight,
        "E": E, "px": px, "py": py, "pz": pz,
        "p": p, "pt": pt, "eta": eta, "phi": phi,
        "gamma": gamma, "beta": beta, "beta_gamma": beta_gamma,
        "mass": np.full_like(p, m_N),
    }


def _empty_dict():
    """Return dict of empty arrays for the no-data case."""
    empty = np.empty(0, dtype=np.float64)
    keys = [
        "weight", "E", "px", "py", "pz",
        "p", "pt", "eta", "phi",
        "gamma", "beta", "beta_gamma", "mass",
    ]
    return {k: empty for k in keys}
