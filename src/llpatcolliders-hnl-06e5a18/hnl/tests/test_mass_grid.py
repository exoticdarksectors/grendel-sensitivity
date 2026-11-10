"""Filename encoding round-trip and uniqueness for the HNL mass grid.

Guards against the historical two-decimal label drift: `mN_0p34` was emitted
for the true grid point 0.335 GeV (5 MeV off). The three-decimal encoding must
be lossless and collision-free across every grid point.
"""

import pytest

from config_mass_grid import (
    MASS_GRID,
    format_mass_for_filename,
    parse_mass_from_filename,
)


@pytest.mark.parametrize("mass", MASS_GRID)
def test_label_roundtrip_is_lossless(mass):
    label = format_mass_for_filename(mass)
    assert parse_mass_from_filename(label) == pytest.approx(mass, abs=1e-9)


@pytest.mark.parametrize("mass", MASS_GRID)
def test_parse_accepts_mN_prefix(mass):
    label = format_mass_for_filename(mass)
    assert parse_mass_from_filename(f"mN_{label}") == pytest.approx(mass, abs=1e-9)


def test_labels_are_unique():
    labels = [format_mass_for_filename(m) for m in MASS_GRID]
    assert len(set(labels)) == len(MASS_GRID)


def test_no_label_drift_exceeds_one_mev():
    # The .2f bug drifted some points by 5 MeV; .3f must keep every label
    # within 1 MeV (one ULP of the encoding) of the true grid value.
    for m in MASS_GRID:
        recovered = parse_mass_from_filename(format_mass_for_filename(m))
        assert abs(recovered - m) <= 0.5e-3
