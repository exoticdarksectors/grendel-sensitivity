"""BC5's production-mode split: the channel column, the per-mode yield curves and their projection."""
import numpy as np
import pandas as pd
import pytest

from grendel.constants import L_INT_PB
from grendel.io.paths import ModelPaths
from grendel.io.variants import split_variant, variants_in
from grendel.io.vectors import csv_column_count, load_combined_csv, write_llp_csv
from grendel.models.scalar import production_bc5 as production
from grendel.models.scalar.spec_bc5 import (ModeYields, QuarticScalarSpec, br_tag, edge_yields,
                                            variant_coefficients)
from grendel.reco.acceptance import scan_u2
from grendel.scan import ScanArrays, ScanConfig, _config_fingerprint


def test_llp_csv_round_trips_the_channel_column(tmp_path):
    n = 5
    rng = np.random.default_rng(0)
    w, E, px, py, pz = rng.random((5, n)) + 1.0
    power = np.array([1, 1, 0, 0, 0], float)
    channel = np.array([0, 0, 1, 2, 2], float)
    path = tmp_path / "mS_1p000.csv"

    write_llp_csv(path, w, E, px, py, pz, coupling_power=power, channel=channel)
    assert csv_column_count(path) == production.N_COLUMNS
    assert production.is_current_format(path)
    data = load_combined_csv(path, 1.0)
    assert np.array_equal(data["coupling_power"], power)
    assert np.array_equal(data["channel"], channel)

    write_llp_csv(path, w, E, px, py, pz, coupling_power=power)
    assert load_combined_csv(path, 1.0)["channel"] is None
    assert not production.is_current_format(path)

    write_llp_csv(path, w, E, px, py, pz, channel=channel)
    data = load_combined_csv(path, 1.0)
    assert np.all(data["coupling_power"] == 1.0) and np.array_equal(data["channel"], channel)

    path.write_text("")
    assert production.is_current_format(path)
    assert load_combined_csv(path, 1.0)["channel"] is None


def test_bc5_rows_carry_their_production_mode():
    rngs = production.make_streams(7)
    pool, higgs, sigma_bottom = production.shared_pools(rngs, 400, 400)
    (w, _, _, _, _, power, channel), counts = production.generate_bc5_4vectors(
        1.0, rngs, pool=pool, higgs=higgs, sigma_bottom=sigma_bottom)
    assert len(w) == len(power) == len(channel) == sum(counts.values())
    assert all(counts[name] > 0 for name in production.CHANNELS)
    for name in production.CHANNELS:
        rows = channel == production.CHANNEL_ID[name]
        assert rows.sum() == counts[name]
        assert np.all(power[rows] == (1.0 if name == "mixing" else 0.0))

    (_, *_, channel_high), counts_high = production.generate_bc5_4vectors(
        20.0, rngs, pool=pool, higgs=higgs, sigma_bottom=sigma_bottom)
    assert counts_high["mixing"] == 0 and counts_high["BSS"] == 0
    assert counts_high["hSS"] == 2 * 400 and np.all(channel_high == production.CHANNEL_ID["hSS"])


def _synthetic_arrays(with_channel=True, n=60, n_samples=4, seed=3):
    rng = np.random.default_rng(seed)
    grid = np.logspace(-12.0, -4.0, 50)
    channel = rng.integers(0, 3, n).astype(float)
    power = np.where(channel == 0, 1.0, 0.0)
    path = rng.uniform(0.5, 3.0, n)
    d = rng.uniform(20.0, 24.0, (n, n_samples))
    passed = rng.random((n, n_samples)) < 0.7
    weights = rng.uniform(1e-3, 1.0, n) * np.where(channel == 0, 1e8, 1.0)
    beta_gamma = rng.uniform(2.0, 50.0, n)
    ctau_ref = 1e-7
    _, N = scan_u2(d, passed, path, weights, beta_gamma, ctau_ref, L_INT_PB, grid,
                   coupling_power=power)
    return ScanArrays(grid=grid, N=N, d=d, passed=passed, path=path, weights=weights,
                      beta_gamma=beta_gamma, ctau_ref=ctau_ref, coupling_power=power,
                      channel=channel if with_channel else None)


def test_every_curve_is_a_linear_combination_of_the_mode_yields():
    arrays = _synthetic_arrays()
    yields = ModeYields(arrays)
    assert set(yields.masks) == {"mixing", "hSS", "BSS"}
    assert np.allclose(sum(yields.on_grid.values()), arrays.N, rtol=1e-12)

    coefficients = variant_coefficients(yields.masks, 0.01, (0.001,))
    assert set(coefficients) == {"mixing", "hSS", "BSS", "brhss0.001", "brhss0.001_hSS",
                                 "brhss0.001_BSS"}
    N, evaluate = yields.combine(coefficients["brhss0.001"])
    assert np.allclose(N, yields.on_grid["mixing"]
                       + 0.1 * (yields.on_grid["hSS"] + yields.on_grid["BSS"]))
    factor = np.where(arrays.channel == 0, 1.0, 0.1)
    _, direct = scan_u2(arrays.d, arrays.passed, arrays.path, arrays.weights * factor,
                        arrays.beta_gamma, arrays.ctau_ref, L_INT_PB, arrays.grid,
                        coupling_power=arrays.coupling_power)
    assert np.allclose(N, direct)

    u2 = 3e-9
    at = yields.at(u2)
    assert np.isclose(evaluate(u2), at["mixing"] + 0.1 * (at["hSS"] + at["BSS"]))
    assert np.isclose(sum(at.values()), arrays.evaluate(u2))
    edges = edge_yields(yields, u2)
    assert edges["mixing"] == at["mixing"]
    assert np.isclose(edges["quartic"], at["hSS"] + at["BSS"])
    assert all(np.isnan(v) for v in edge_yields(yields, np.nan).values())


def test_without_the_channel_column_only_mixing_and_quartic_are_told_apart():
    arrays = _synthetic_arrays(with_channel=False)
    yields = ModeYields(arrays)
    assert set(yields.masks) == {"mixing", "quartic"}
    coefficients = variant_coefficients(yields.masks, 0.01, (0.001, 0.1))
    assert set(coefficients) == {"mixing", "quartic", "brhss0.001", "brhss0.001_quartic",
                                 "brhss0.1", "brhss0.1_quartic"}
    assert coefficients["brhss0.1"] == {"mixing": 1.0, "quartic": 10.0}
    edges = edge_yields(yields, 1e-9)
    assert np.isnan(edges["hSS"]) and np.isnan(edges["BSS"]) and edges["quartic"] > 0.0


def test_br_tag():
    assert br_tag(0.001) == "brhss0.001"
    assert br_tag(0.01) == "brhss0.01"
    assert br_tag(1e-5) == "brhss1e-05"


def test_spec_validates_the_overlays_and_records_them(tmp_path):
    paths = ModelPaths.resolve("bc5", vectors=tmp_path, templates=tmp_path, analysis=tmp_path,
                               geometry=tmp_path)
    spec = QuarticScalarSpec(paths, br_hss_overlays=(0.001, 0.1))
    assert spec.br_hss_overlays == (0.001, 0.1)
    assert QuarticScalarSpec(paths, br_hss=0.0).br_hss_overlays == ()
    with pytest.raises(ValueError):
        QuarticScalarSpec(paths, br_hss_overlays=(0.0,))
    fingerprint = _config_fingerprint(spec, ScanConfig(decay_samples=1))
    assert fingerprint["model"]["br_hss_overlays"] == [0.001, 0.1]


def test_split_variant_puts_one_variant_in_the_nominal_place():
    frame = pd.DataFrame({
        "mass_GeV": [1.0, 2.0], "n_hits": [10, 20],
        "u2_min": [1e-9, 2e-9], "u2_max": [1e-6, 2e-6], "peak_N": [50.0, 60.0],
        "has_sensitivity": [True, True], "u2_min_N10": [2e-9, 3e-9],
        "u2_min_N_mixing": [1.0, 2.0], "u2_min_sample_ess": [5.0, 6.0],
        "mixing_u2_min": [3e-9, np.nan], "mixing_u2_max": [5e-7, np.nan],
        "mixing_peak_N": [20.0, 1.0], "mixing_has_sensitivity": [True, False],
        "mixing_u2_min_N10": [4e-9, np.nan],
        "hSS_u2_min": [4e-9, 2.5e-9], "hSS_u2_max": [8e-7, 1.5e-6], "hSS_peak_N": [30.0, 59.0],
        "hSS_has_sensitivity": [True, True], "hSS_u2_min_N10": [5e-9, 3.5e-9],
        "brhss0.001_hSS_u2_min": [6e-9, 4e-9],
    })
    assert variants_in(frame) == ["mixing", "hSS", "brhss0.001_hSS"]
    out = split_variant(frame, "hSS")
    assert list(out.columns) == ["mass_GeV", "n_hits", "u2_min", "u2_max", "peak_N",
                                 "has_sensitivity", "u2_min_N10", "u2_min_N_mixing",
                                 "u2_min_sample_ess"]
    assert out["u2_min"].tolist() == frame["hSS_u2_min"].tolist()
    assert out["u2_min_N10"].tolist() == frame["hSS_u2_min_N10"].tolist()
    assert out["n_hits"].tolist() == frame["n_hits"].tolist()
    with pytest.raises(ValueError, match="variants present: mixing, hSS"):
        split_variant(frame, "BSS")


def test_scan_u2_scans_a_weight_matrix_row_by_row():
    arrays = _synthetic_arrays()
    masks = [arrays.channel == i for i in range(3)]
    matrix = np.stack([arrays.weights * m for m in masks])
    for power in (arrays.coupling_power, None):
        _, together = scan_u2(arrays.d, arrays.passed, arrays.path, matrix, arrays.beta_gamma,
                              arrays.ctau_ref, L_INT_PB, arrays.grid, coupling_power=power)
        assert together.shape == (3, len(arrays.grid))
        for row, m in zip(together, masks):
            _, alone = scan_u2(arrays.d, arrays.passed, arrays.path, arrays.weights * m,
                               arrays.beta_gamma, arrays.ctau_ref, L_INT_PB, arrays.grid,
                               coupling_power=power)
            assert alone.shape == (len(arrays.grid),) and np.allclose(row, alone, rtol=1e-13)
