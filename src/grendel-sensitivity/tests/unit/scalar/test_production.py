import numpy as np
import pytest

from grendel.models.scalar import model
from grendel.models.scalar.production import B_SPECIES, generate_scalar_4vectors


def test_scalar_production_propagates_fonll_importance_weights():
    sampling_weight = np.array([0.5, 1.0, 1.5])
    pool = {
        "pt": np.array([2.0, 4.0, 8.0]),
        "y": np.array([-0.2, 0.0, 0.3]),
        "phi": np.array([0.1, 1.2, 2.4]),
        "sampling_weight": sampling_weight,
    }

    weights, *_ = generate_scalar_4vectors(
        0.5,
        len(sampling_weight),
        np.random.default_rng(7),
        sigma_bottom=100.0,
        pool=pool,
    )

    by_species = weights.reshape(len(B_SPECIES), len(sampling_weight))
    expected_shape = sampling_weight / sampling_weight.mean()
    for row in by_species:
        assert np.allclose(row / row.mean(), expected_shape)


def test_full_b_hadron_pool_including_baryons():
    # The four species (B+, B0, Bs, Lambda_b) must span the whole b-hadron pool.
    tags = {tag for (tag, *_rest) in B_SPECIES.values()}
    assert tags == {"B+", "B0", "Bs", "Lambda_b"}
    frag_sum = sum(frag for (*_h, frag) in B_SPECIES.values())
    assert frag_sum == pytest.approx(1.0, abs=1e-6)


def test_lambda_b_uses_the_lambda_recoil_not_the_kaon():
    # The recoil mass sets the kinematic ceiling. With the physical Lambda recoil
    # Lambda_b closes at m_Lambda_b - m_Lambda ~ 4.50 GeV; a kaon recoil would
    # wrongly leak it to ~5.13 GeV. Check the entry carries the Lambda mass and
    # that production is closed for Lambda_b between the two ceilings.
    _tag, m_parent, m_recoil, _frag = B_SPECIES[5122]
    assert m_recoil == pytest.approx(model.M_LAMBDA, abs=1e-6)
    assert m_parent - m_recoil == pytest.approx(4.50, abs=0.02)

    # At m_S = 4.7 GeV only the mesons (ceilings ~4.78-4.87) remain; Lambda_b is
    # closed. Its BR is still positive there, so a kaon recoil would have kept it.
    rng = np.random.default_rng(0)
    pool = {"pt": np.array([3.0]), "y": np.array([0.0]), "phi": np.array([0.0])}
    assert float(model.br_B_to_Xs_S(4.7, parent="Lambda_b")) > 0.0
    w, *_ = generate_scalar_4vectors(4.7, 1, rng, sigma_bottom=1.0, pool=pool)
    # exactly the meson species (3) contribute one event each at this mass
    assert len(w) == 3


def test_baryons_raise_the_low_mass_yield_by_the_expected_fraction():
    # Adding the b-baryon pool is an additive production boost; at the deep-reach
    # mass (0.975 GeV) it is recoil-insensitive and should lift the total weight
    # by frag_Lambda_b * BR_Lambda_b / sum_meson(frag * BR) ~ +23%.
    m_S = 0.975
    rng = np.random.default_rng(3)
    pool = {"pt": np.array([2.0, 5.0]), "y": np.array([0.0, 0.5]),
            "phi": np.array([0.3, 1.1])}
    total = 0.0
    baryon = 0.0
    for pdg, (parent, m_B, m_recoil, frag) in B_SPECIES.items():
        br = float(model.br_B_to_Xs_S(m_S, parent=parent))
        contrib = frag * br
        total += contrib
        if pdg == 5122:
            baryon = contrib
    meson = total - baryon
    assert baryon / meson == pytest.approx(0.23, abs=0.03)
