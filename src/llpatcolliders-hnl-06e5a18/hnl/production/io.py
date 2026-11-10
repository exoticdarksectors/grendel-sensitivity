"""CSV and path helpers shared by HNL production channels."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from config_mass_grid import format_mass_for_filename
from production.paths import LLP_VECTORS_DIR

CSV_FMT = "%.8e"


def llp_csv_path(flavor, channel, mass, *, base=LLP_VECTORS_DIR, mkdir=True) -> Path:
    path = base / flavor / channel / f"mN_{format_mass_for_filename(mass)}.csv"
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
