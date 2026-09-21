"""Parse FONLL meson-level differential cross-section tables.

The bottom and charm grids are selected independently so a single production
run can use, e.g., a varied m_b grid for bottom while holding charm central.
Each quark's grid is overridden by its own environment variable, pointing at
a *coherent* individual grid (one consistent scale/PDF-member/mass choice) --
never a pointwise envelope grid, which is not a physical cross section:

    GRENDEL_FONLL_BOTTOM_GRID   path to the bottom dsigma/dpT/dy table
    GRENDEL_FONLL_CHARM_GRID    path to the charm  dsigma/dpT/dy table

Unset variables fall back to the central grids shipped in ``data/fonll/central``
(see the PROVENANCE file there for how they were generated).
"""

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
    """Resolve the FONLL grid for ``quark``: per-quark env override, else central.

    The override must exist; a missing path is an error rather than a silent
    fall-back, so a mistyped variation grid cannot quietly revert to central.
    """
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
    """Refuse pointwise variation-envelope grids.

    The variation combiner (``grendel.production.fonll_grids.combine``)
    writes envelope grids (scale/pdf/mass max-min or +/-sigma constructions)
    in this same three-column format and
    marks them with an ``envelope_band`` header. They are not physical cross
    sections, so sampling kinematics from them is meaningless; only coherent
    individual variation grids may be used here.
    """
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
    integral_over_y = trapezoid(dsigma_2d, y_arr, axis=1)  # shape (N_pt,)
    sigma_total = trapezoid(integral_over_y, pt_arr)        # scalar (pb)

    return sigma_total
