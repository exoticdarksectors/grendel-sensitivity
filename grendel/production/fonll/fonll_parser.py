"""Parse FONLL meson-level differential cross-section tables."""

import os
import numpy as np
from pathlib import Path

_CENTRAL_DIR = Path(__file__).resolve().parent.parent / "data" / "fonll" / "central"

_DEFAULT_GRIDS = {
    "bottom": _CENTRAL_DIR / "fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_bottom.dat",
    "charm": _CENTRAL_DIR / "fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_charm.dat",
}

_GRID_ENV_VARS = {
    "bottom": "GRENDEL_FONLL_BOTTOM_GRID",
    "charm": "GRENDEL_FONLL_CHARM_GRID",
}


def fonll_grid_path(quark):
    """Resolve the FONLL grid for ``quark``: per-quark env override, else central."""
    if quark not in _DEFAULT_GRIDS:
        raise KeyError(f"unknown quark {quark!r}; expected one of {sorted(_DEFAULT_GRIDS)}")
    override = os.environ.get(_GRID_ENV_VARS[quark])
    if override:
        path = Path(override).expanduser()
        if not path.exists():
            raise FileNotFoundError(
                f"{_GRID_ENV_VARS[quark]}={override} does not exist"
            )
        return path
    return _DEFAULT_GRIDS[quark]


class _GridFileMap:
    """Backward-compatible ``FONLL_FILES[quark]`` that resolves at access time."""

    def __getitem__(self, quark):
        return fonll_grid_path(quark)

    def __contains__(self, quark):
        return quark in _DEFAULT_GRIDS

    def keys(self):
        return _DEFAULT_GRIDS.keys()

    def items(self):
        return [(quark, fonll_grid_path(quark)) for quark in _DEFAULT_GRIDS]


FONLL_FILES = _GridFileMap()


def _reject_envelope_grid(path):
    """Refuse pointwise variation-envelope grids."""
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if "envelope_band:" in line:
                raise ValueError(
                    f"{path} is a pointwise variation envelope, not a physical "
                    f"cross section; point the grid override at a coherent "
                    f"individual variation grid (one scale point, one PDF "
                    f"member, one mass) instead"
                )


def parse_fonll_file(path):
    _reject_envelope_grid(path)
    data = np.loadtxt(path, comments="#")
    pt_all = data[:, 0]
    y_all = data[:, 1]
    dsigma_all = data[:, 2]

    pt_unique = np.unique(pt_all)
    y_unique = np.unique(y_all)
    n_pt = len(pt_unique)
    n_y = len(y_unique)

    if len(data) != n_pt * n_y:
        raise ValueError(f"{path}: expected a rectangular pT-y grid")
    expected_pt = np.repeat(pt_unique, n_y)
    expected_y = np.tile(y_unique, n_pt)
    if not np.array_equal(pt_all, expected_pt) or not np.array_equal(y_all, expected_y):
        if np.array_equal(y_all[:n_y], y_unique[::-1]):
            raise ValueError(f"{path}: y column must be in ascending order")
        raise ValueError(f"{path}: rows must be ordered with pT outermost and y innermost")
    dsigma_2d = dsigma_all.reshape(n_pt, n_y)

    return pt_unique, y_unique, dsigma_2d


def get_sigma_total(quark):
    path = fonll_grid_path(quark)
    pt_arr, y_arr, dsigma_2d = parse_fonll_file(path)

    trapezoid = getattr(np, "trapezoid", None) or np.trapz
    integral_over_y = trapezoid(dsigma_2d, y_arr, axis=1)
    sigma_total = trapezoid(integral_over_y, pt_arr)

    return sigma_total
