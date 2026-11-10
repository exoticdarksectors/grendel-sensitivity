"""Sanity checks on the FONLL table parser using the vendored FONLL tables."""

import numpy as np
import pytest

from production.fonll.fonll_parser import (
    FONLL_FILES, parse_fonll_file, get_sigma_total,
)


def test_table_files_exist():
    for q, path in FONLL_FILES.items():
        assert path.exists(), f"Missing vendored FONLL table for {q}: {path}"


def test_parse_shape_is_rectangular():
    """Tables parse as 1D pT, 1D y, 2D dsigma of consistent shape.
    Does NOT hard-code grid density."""
    for q in ("charm", "bottom"):
        pt, y, ds = parse_fonll_file(FONLL_FILES[q])
        assert pt.ndim == 1
        assert y.ndim == 1
        assert ds.shape == (len(pt), len(y))
        assert len(pt) >= 2 and len(y) >= 2
        assert pt.min() <= 0.0 and pt.max() >= 50.0
        assert y.min() <= -3.0 and y.max() >= 3.0


def test_nnpdf40_grid_matches_dense_production_contract():
    for q in ("charm", "bottom"):
        pt, y, ds = parse_fonll_file(FONLL_FILES[q])
        assert pt.ndim == 1
        assert y.ndim == 1
        assert ds.shape == (100, 100)
        assert np.allclose(pt, np.linspace(0.0, 50.0, 100))
        assert np.allclose(y, np.linspace(-3.0, 3.0, 100))
        assert np.isfinite(ds).all()


def test_parser_rejects_rows_outside_pt_outer_y_inner_order(tmp_path):
    path = tmp_path / "shuffled.dat"
    np.savetxt(path, np.array([
        [0.0, -1.0, 1.0],
        [1.0, -1.0, 2.0],
        [0.0, 1.0, 3.0],
        [1.0, 1.0, 4.0],
    ]))
    with pytest.raises(ValueError, match="pT outermost and y innermost"):
        parse_fonll_file(path)


def test_parser_rejects_descending_y(tmp_path):
    """A file with descending y rows must raise a clear ascending-order error."""
    path = tmp_path / "descending_y.dat"
    pt_nodes = np.array([0.0, 1.0])
    y_nodes_desc = np.array([1.0, 0.0, -1.0])  # descending
    rows = []
    for pt in pt_nodes:
        for y in y_nodes_desc:
            rows.append([pt, y, 1.0])
    np.savetxt(path, np.array(rows))
    with pytest.raises(ValueError, match="y column must be in ascending order"):
        parse_fonll_file(path)


def test_parser_rejects_envelope_grids(tmp_path):
    """Pointwise variation envelopes are not physical cross sections and must
    be refused even though they share the three-column format."""
    path = tmp_path / "envelope.dat"
    rows = []
    for pt in (0.0, 1.0):
        for y in (-1.0, 0.0, 1.0):
            rows.append([pt, y, 1.0])
    header = "# FONLL heavy-flavor meson grid (variation envelope)\n# envelope_band: scale_up\n"
    path.write_text(header + "\n".join(" ".join(str(v) for v in r) for r in rows) + "\n")
    with pytest.raises(ValueError, match="variation envelope"):
        parse_fonll_file(path)


def test_dsigma_nonnegative_majority():
    """At least 95% of bins should be non-negative (small numerical negatives allowed near edges)."""
    for q in ("charm", "bottom"):
        _, _, ds = parse_fonll_file(FONLL_FILES[q])
        frac_nn = np.mean(ds >= 0)
        assert frac_nn > 0.95, f"{q}: only {frac_nn:.1%} non-negative bins"


def test_total_sigma_order_of_magnitude():
    """Total integrated cross section should match expected LHC 14 TeV ballpark."""
    # FONLL at 14 TeV: sigma_cc ~ 1e10 pb, sigma_bb ~ 1e8 pb (very approximate)
    sigma_c = get_sigma_total("charm")
    sigma_b = get_sigma_total("bottom")
    assert 1e9 < sigma_c < 1e11, f"sigma_charm = {sigma_c:.3e} pb out of expected range"
    assert 1e7 < sigma_b < 1e9, f"sigma_bottom = {sigma_b:.3e} pb out of expected range"
    # charm > bottom at LHC
    assert sigma_c > sigma_b
