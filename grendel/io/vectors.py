"""Four-vector files and their filename convention.

A production channel writes one headerless CSV per mass point with columns
``weight, E, px, py, pz`` (weight in pb at unit coupling, momenta in GeV).
``load_combined_csv`` turns such a file into the derived kinematics the
ray-cast and the acceptance Monte Carlo consume.

Mass labels carry three decimals with the decimal point written as ``p``
(``1.025 -> 1p025``): the mass grids have 15-MeV and 25-MeV spacings, so a
two-decimal encoding would alias neighbouring points.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

CSV_FMT = "%.8e"


def format_mass_for_filename(mass) -> str:
    """Filename encoding: 1.025 -> '1p025'."""
    return f"{mass:.3f}".replace(".", "p")


def parse_mass_from_filename(label: str) -> float:
    """Inverse of ``format_mass_for_filename``; accepts a bare label
    (``'1p025'``) or a stem with a one-letter prefix (``'mN_1p025'``)."""
    stem = label.split("_", 1)[1] if "_" in label else label
    return float(stem.replace("p", "."))


def llp_csv_path(base, flavor, channel, mass, *, prefix="mN", mkdir=True) -> Path:
    """``<base>/<flavor>/<channel>/<prefix>_<label>.csv`` for one mass point."""
    path = Path(base) / flavor / channel / f"{prefix}_{format_mass_for_filename(mass)}.csv"
    if mkdir:
        path.parent.mkdir(parents=True, exist_ok=True)
    return path


def read_csv_matrix(path) -> np.ndarray:
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return np.empty((0, 0))
    data = np.loadtxt(path, delimiter=",")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def write_csv_matrix(path, data) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if len(data) == 0:
        path.write_text("")
    else:
        np.savetxt(path, data, delimiter=",", fmt=CSV_FMT)


def write_llp_csv(path, weights, E, px, py, pz) -> None:
    write_csv_matrix(path, np.column_stack([weights, E, px, py, pz]))


def write_empty_csv(path) -> None:
    write_csv_matrix(path, [])


def scale_weight_column(path, factor) -> None:
    data = read_csv_matrix(path)
    if data.size:
        data[:, 0] *= factor
        write_csv_matrix(path, data)


def load_combined_csv(csv_path, mass):
    """Load a four-vector CSV and compute the derived kinematics.

    Returns a dict of equal-length arrays: ``weight``, ``E``, ``px``, ``py``,
    ``pz``, ``p``, ``pt``, ``eta``, ``phi``, ``gamma``, ``beta``,
    ``beta_gamma`` and ``mass`` (filled with ``mass``). Empty arrays when the
    file is missing or empty.
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

    theta = np.arctan2(pt, pz)
    theta = np.clip(theta, 1e-10, np.pi - 1e-10)
    eta = -np.log(np.tan(theta / 2.0))

    phi = np.arctan2(py, px)

    energy = np.maximum(E, mass)  # guard against numerical noise
    gamma = energy / mass
    beta = p / energy
    beta_gamma = p / mass

    return {
        "weight": weight,
        "E": E, "px": px, "py": py, "pz": pz,
        "p": p, "pt": pt, "eta": eta, "phi": phi,
        "gamma": gamma, "beta": beta, "beta_gamma": beta_gamma,
        "mass": np.full_like(p, mass),
    }


def _empty_dict():
    empty = np.empty(0, dtype=np.float64)
    keys = ["weight", "E", "px", "py", "pz", "p", "pt", "eta", "phi",
            "gamma", "beta", "beta_gamma", "mass"]
    return {k: empty for k in keys}
