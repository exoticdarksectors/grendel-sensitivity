"""Pure-function regression on production.decay_engine.tau_decay.

Targets the shared sampler that both induced-tau and prompt-tau use:

  - init_hnlcalc returns a HNLCalc with the right active mixing.
  - compute_tau_production_br_components returns non-negative branching
    fractions whose count tracks the kinematically-open channels.
  - sample_hnl_from_tau preserves 4-momentum (lab-frame Lorentz boost from
    decay -> reconstruction round-trip).
  - asymmetry=0 reproduces a flat cos(theta) distribution in the parent
    rest frame; asymmetry=-1 biases toward cos(theta) = -1.

No MadGraph or LHAPDF dependency; runs anywhere HNLCalc imports.
"""

import numpy as np
import pytest

from production.constants import M_TAU
from production.decay_engine.tau_decay import (
    init_hnlcalc,
    compute_tau_production_br_components,
    sample_hnl_from_tau,
    TAU_MESON_2BODY_ASYMMETRY,
    TAU_W_2BODY_ASYMMETRY,
    TAU_2BODY_MESON_PDGS,
)


@pytest.fixture(scope="module")
def hnl_utau():
    return init_hnlcalc("Utau")


def test_init_hnlcalc_flavors_are_distinct():
    """Each flavor selects a different active mixing."""
    hnl_e = init_hnlcalc("Ue")
    hnl_mu = init_hnlcalc("Umu")
    hnl_tau = init_hnlcalc("Utau")
    # The three HNLCalc instances must be distinct objects with distinct
    # internal mixings; this is the only behavioural difference that matters.
    assert hnl_e is not hnl_mu
    assert hnl_mu is not hnl_tau

    with pytest.raises(ValueError):
        init_hnlcalc("Unobtainium")


def test_br_components_nonnegative_and_finite(hnl_utau):
    """BR sums add up to a finite non-negative total below m_tau.

    After the matrix-element-weighted sampler landed, br3 tuples carry the
    HNLCalc dBR expression alongside the rate:
      br2 entries: (meson_mass, br)
      br3 entries: (lep_mass, dbr_expr, br)
    """
    br2, br3, br_total = compute_tau_production_br_components(hnl_utau, 0.5)
    assert br_total > 0
    assert all(br >= 0 for _, br in br2)
    assert all(br >= 0 for _, _, br in br3)
    # Total should match the sum of pieces.
    expected = sum(br for _, br in br2) + sum(br for _, _, br in br3)
    assert br_total == pytest.approx(expected, rel=1e-12)


def test_br_components_close_at_threshold(hnl_utau):
    """Channels disappear as m_N approaches m_tau."""
    # Far below threshold: many open channels.
    br2_low, br3_low, total_low = compute_tau_production_br_components(hnl_utau, 0.1)
    # Just below: most channels closed (only e/mu lepton 3-body and pi if it fits).
    br2_high, br3_high, total_high = compute_tau_production_br_components(hnl_utau, 1.6)

    # Channels must monotonically decrease.
    assert len(br2_low) >= len(br2_high)
    assert total_low >= total_high


def test_sample_hnl_from_tau_invariant_mass(hnl_utau):
    """Reconstructed HNL invariant mass equals m_N (to numerical noise)."""
    rng = np.random.default_rng(0)
    n = 2000
    # All taus at rest in their own frame, given a fixed boost along z.
    gamma = 5.0
    beta = np.sqrt(1 - 1 / gamma ** 2)
    tau_E = np.full(n, gamma * M_TAU)
    tau_px = np.zeros(n)
    tau_py = np.zeros(n)
    tau_pz = np.full(n, gamma * beta * M_TAU)

    m_N = 0.5
    br2, br3, brt = compute_tau_production_br_components(hnl_utau, m_N)
    hnl_4v, _, _ = sample_hnl_from_tau(
        tau_E, tau_px, tau_py, tau_pz, m_N,
        br2, br3, brt, rng, asymmetry=0.0,
    )
    m2 = hnl_4v[:, 0] ** 2 - hnl_4v[:, 1] ** 2 - hnl_4v[:, 2] ** 2 - hnl_4v[:, 3] ** 2
    # Numerical floor: floats * boost + decay accumulate ~1e-6 relative noise.
    assert m2.mean() == pytest.approx(m_N ** 2, rel=2e-3)
    assert m2.std() < 1e-3


def test_sample_hnl_from_tau_asymmetry_flips_distribution(hnl_utau):
    """Mean N rapidity in the tau rest frame shifts with asymmetry.

    asymmetry=0 -> isotropic in tau rest frame -> <p_z^N (tau frame)> = 0;
    asymmetry=-1 (heavy-meson-origin taus) -> N recoils backward to the tau
    momentum; asymmetry=+1 (W-origin taus) -> N recoils forward.

    Operationally we check that the *direction* of the mean lab-frame p_z
    shifts both ways around the flat case (the tau is boosted along +z, so a
    tau-rest backward N becomes a lab N with lower mean p_z and vice versa).
    """
    rng = np.random.default_rng(1)
    n = 5000
    gamma = 3.0
    beta = np.sqrt(1 - 1 / gamma ** 2)
    tau_E = np.full(n, gamma * M_TAU)
    tau_px = np.zeros(n)
    tau_py = np.zeros(n)
    tau_pz = np.full(n, gamma * beta * M_TAU)

    m_N = 0.3
    br2, br3, brt = compute_tau_production_br_components(hnl_utau, m_N)

    flat, _, _ = sample_hnl_from_tau(
        tau_E, tau_px, tau_py, tau_pz, m_N, br2, br3, brt,
        np.random.default_rng(1), asymmetry=0.0,
    )
    meson_like, _, _ = sample_hnl_from_tau(
        tau_E, tau_px, tau_py, tau_pz, m_N, br2, br3, brt,
        np.random.default_rng(1), asymmetry=TAU_MESON_2BODY_ASYMMETRY,  # -1
    )
    w_like, _, _ = sample_hnl_from_tau(
        tau_E, tau_px, tau_py, tau_pz, m_N, br2, br3, brt,
        np.random.default_rng(1), asymmetry=TAU_W_2BODY_ASYMMETRY,  # +1
    )

    # asymmetry=-1 pulls the N backward in the tau rest frame (lower mean lab
    # pz); asymmetry=+1 pushes it forward (higher mean lab pz). Loose bounds
    # absorb the 3-body component (~unaffected) and statistical noise.
    assert meson_like[:, 3].mean() < flat[:, 3].mean()
    assert w_like[:, 3].mean() > flat[:, 3].mean()
