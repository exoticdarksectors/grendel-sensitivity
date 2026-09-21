"""Unit tests for the exHad decay-template packer (no exHad runtime needed)."""
import numpy as np
import pytest

from grendel.models.hnl.templates import exhad as hx
from grendel.reco import templates as gx


def test_flatten_events_packs_bundle_and_counts_charged():
    # two decays: N -> mu- pi+ (2 charged) and N -> 3 nu (0 charged)
    events = [
        [[0.1, 0.2, 0.3, 0.5, 0.10566, 13.0], [-0.1, -0.2, -0.3, 0.5, 0.13957, 211.0]],
        [[0.0, 0.0, 0.4, 0.4, 0.0, 14.0], [0.0, 0.3, -0.2, 0.36, 0.0, -14.0],
         [0.0, -0.3, -0.2, 0.36, 0.0, 12.0]],
    ]
    b = gx.flatten_events(events, ["Pimu", "3nu"])
    assert b["daughter_counts"].tolist() == [2, 3]
    assert b["pdg"].tolist() == [13, 211, 14, -14, 12]
    assert b["charge"].tolist() == [-1.0, 1.0, 0.0, 0.0, 0.0]
    assert b["stable"].all() and b["stable"].dtype == bool
    assert b["n_charged"].tolist() == [2, 0]
    assert b["channel_label"].tolist() == ["Pimu", "3nu"]
    assert np.allclose(b["energy"], [0.5, 0.5, 0.4, 0.36, 0.36])
    assert np.allclose(b["mass"], [0.10566, 0.13957, 0.0, 0.0, 0.0])


def test_charge_lookup_covers_the_final_state_species_exhad_returns():
    expected = {211: 1, -211: -1, 321: 1, -321: -1, 11: -1, -11: 1, 13: -1, -13: 1,
                2212: 1, -2212: -1, 2112: 0, 130: 0, 22: 0, 12: 0, 14: 0, 16: 0}
    for pdg, q in expected.items():
        assert gx.charge_of(pdg) == q, pdg


def test_ctau_interpolation_is_log_log_and_flavour_specific():
    hbarc = gx.HBARC_GEV_M
    table = np.array([[1.0, 1e-15, 2e-15, 4e-15],
                      [4.0, 1e-13, 2e-13, 4e-13]])   # Gamma ~ m^{10/3}: exact in log-log
    assert hx.ctau_u2eq1_exhad(1.0, "Ue", table) == pytest.approx(hbarc / 1e-15)
    assert hx.ctau_u2eq1_exhad(4.0, "Umu", table) == pytest.approx(hbarc / 2e-13)
    mid = hx.ctau_u2eq1_exhad(2.0, "Utau", table)
    assert mid == pytest.approx(hbarc / (4e-15 * 2.0 ** (np.log(100) / np.log(4))), rel=1e-9)
    with pytest.raises(ValueError):
        hx.ctau_u2eq1_exhad(0.5, "Ue", table)


def test_point_seed_is_deterministic_distinct_and_in_exhad_range():
    a = gx.point_seed(1234, "Ue", "1p000", namespace=hx.SEED_NAMESPACE)
    assert a == gx.point_seed(1234, "Ue", "1p000", namespace=hx.SEED_NAMESPACE)
    assert a != gx.point_seed(1234, "Umu", "1p000", namespace=hx.SEED_NAMESPACE)
    assert a != gx.point_seed(1234, "Ue", "1p025", namespace=hx.SEED_NAMESPACE)
    assert a != gx.point_seed(1235, "Ue", "1p000", namespace=hx.SEED_NAMESPACE)
    assert 0 <= a < 2 ** 63
