"""The exHad decay-model supplement to the BC4 variation envelope."""
import numpy as np
import pytest
import pandas as pd

from grendel.band.scalar import exhad_variation as xv
from grendel.band.campaign_store import stable_seed
from grendel.band.scalar import variations as ub


def test_template_seed_is_derived_from_the_variation_and_distinct():
    seed = xv.templates_seed(42)
    assert seed == xv.templates_seed(42)
    assert seed != stable_seed(42, ub.EXHAD_DECAY_VARIATION, "parent_pool")
    assert seed != stable_seed(42, ub.DECAY_VARIATION, "parent_pool")
    assert 0 <= seed < 2**32


def test_containment_counts_breaches_in_dex():
    band = pd.DataFrame({
        "mass_GeV": [1.0, 2.0, 3.0], "has_sensitivity": [True, True, True],
        "u2_min_envelope_lo": [1e-9, 1e-9, np.nan], "u2_min_envelope_hi": [1e-8, 1e-8, np.nan],
        "u2_max_envelope_lo": [1e-5, 1e-5, 1e-5], "u2_max_envelope_hi": [1e-4, 1e-4, 1e-4],
    })
    curve = pd.DataFrame({
        "mass_GeV": [1.0, 2.0, 3.0], "has_sensitivity": [True, True, False],
        "u2_min": [5e-9, 1e-10, 5e-9], "u2_max": [5e-5, 1e-3, 5e-5],
    })
    out = xv.containment(curve, band)
    assert (out["u2_min"]["n"], out["u2_min"]["inside"]) == (2, 1)
    assert out["u2_min"]["worst_breach_dex"] == pytest.approx(1.0)
    assert out["u2_max"]["n"] == 2 and out["u2_max"]["inside"] == 1
    assert out["u2_max"]["worst_breach_dex"] == pytest.approx(1.0)
