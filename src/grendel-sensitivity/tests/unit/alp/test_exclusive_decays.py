"""Checks for the decoded 2501 exclusive-decay interface."""

import sys
from pathlib import Path

import pytest

from grendel.models.alp import exclusive_decays as decays
from grendel.models.alp import mass_grid
from grendel.models.alp import model


def test_mass_grid_is_shared_dense_scan():
    assert len(mass_grid.ALP_MASS_GRID) == 145
    assert mass_grid.ALP_MASS_GRID[0] == pytest.approx(0.22)
    assert mass_grid.ALP_MASS_GRID[-1] == pytest.approx(4.75)
    assert all(
        mass in mass_grid.ALP_MASS_GRID
        for mass in (1.18, 1.19, 1.60, 3.12, 3.34, 3.38)
    )


@pytest.mark.parametrize("mass", [0.12, 0.22, 1.0, 2.0, 2.5, 4.75])
def test_exclusive_weights_are_normalized_and_unique(mass):
    weights = decays.exclusive_branching_weights(mass)
    assert sum(weights.values()) == pytest.approx(1.0, abs=1e-12)
    assert not decays.DUPLICATE_CHANNEL_IDS.intersection(weights)
    assert all(weight > 0.0 for weight in weights.values())


def test_every_unique_table_channel_has_a_pythia_mapping():
    table_channels = set(decays._branching_table().dtype.names) - {"mass_GeV"}
    expected = table_channels - decays.DUPLICATE_CHANNEL_IDS
    assert expected == set(decays.DECAY_PRODUCTS_PDG) - decays.DUPLICATE_CHANNEL_IDS


@pytest.mark.parametrize("mass", [0.22, 1.0, 2.5, 4.0])
def test_reference_ctau_matches_model(mass):
    assert decays.ctau_at_reference_coupling(mass) == pytest.approx(
        model.alp_ctau(mass, model.INV_F_REF), rel=1e-12
    )
