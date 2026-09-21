"""Template decay backend inputs and the secondary-threshold columns of the
BC4 scan."""
import numpy as np
import pytest

from grendel.models.hnl.spec import load_template_bundle
from grendel.models.scalar.spec import BAND_FIELDS
from grendel.reco.exclusion import find_exclusion_band_refined
from grendel.scan import secondary_threshold_columns


def test_template_bundle_loads_by_mass_label(tmp_path):
    np.savez_compressed(tmp_path / "templates_0p975.npz",
                        daughter_counts=np.array([2], np.int32), pdg=np.array([211, -211], np.int32),
                        px=np.zeros(2), py=np.zeros(2), pz=np.array([0.4, -0.4]),
                        energy=np.array([0.4875, 0.4875]), mass=np.array([0.13957, 0.13957]),
                        charge=np.array([1.0, -1.0]), stable=np.ones(2, bool),
                        mass_GeV=np.float64(0.975), ctau_m_u2eq1=np.float64(2.6e-10),
                        decay_model=np.array("exhad:scalar-1809"))
    bundle = load_template_bundle(tmp_path / "templates_0p975.npz")
    assert bundle is not None
    assert float(bundle["ctau_m_u2eq1"]) == pytest.approx(2.6e-10)
    assert bundle["daughter_counts"].tolist() == [2]
    assert str(bundle["decay_model"]) == "exhad:scalar-1809"
    assert load_template_bundle(tmp_path / "templates_1p000.npz") is None


def test_secondary_thresholds_add_consistent_columns():
    # A lognormal-like yield curve peaking at 1e-9 with peak 40: the N>=3
    # island is wider than the N>=10 one and both are closed.
    grid = np.logspace(-12, -2, 200)
    N = 40.0 * np.exp(-0.5 * (np.log10(grid) + 9.0) ** 2 / 0.6 ** 2)

    def evaluate(u2):
        return float(np.interp(np.log10(u2), np.log10(grid), N))

    def solve(threshold):
        return find_exclusion_band_refined(grid, N, evaluate, threshold)

    result = {}
    secondary_threshold_columns(result, solve, (10.0,), BAND_FIELDS)
    assert result["has_sensitivity_N10"]
    assert result["u2_min_N10"] > 0 and result["u2_max_N10"] > result["u2_min_N10"]
    assert not result["u2_min_open_N10"] and not result["u2_max_open_N10"]
    assert evaluate(result["u2_min_N10"]) == pytest.approx(10.0, rel=1e-3)
    assert evaluate(result["u2_max_N10"]) == pytest.approx(10.0, rel=1e-3)
    base = solve(3.0)
    assert base["u2_min"] < result["u2_min_N10"] < result["u2_max_N10"] < base["u2_max"]
    empty = {}
    secondary_threshold_columns(empty, solve, (), BAND_FIELDS)
    assert empty == {}
