"""Sampled meson 4-vectors should be on-shell and respect the FONLL grid bounds."""

import numpy as np
import pytest

from production.constants import MESON_MASSES, FRAG_B, FRAG_C, OMITTED_FRAG_B, OMITTED_FRAG_C
from production.decay_engine.generate_meson_csvs import generate_bc_pool
from production.fonll.meson_sampler import (
    _build_cdf, _node_bin_edges, _sample_node_intervals, sample_meson_4vectors,
)


def test_edge_bin_smearing_stays_inside_table_without_endpoint_clipping():
    rng = np.random.default_rng(0)
    edges = _node_bin_edges(np.array([0.0, 1.0, 2.0]))
    indices = np.repeat([0, 2], 1000)
    sampled = _sample_node_intervals(edges, indices, rng)

    assert np.all((sampled[:1000] >= 0.0) & (sampled[:1000] <= 0.5))
    assert np.all((sampled[1000:] >= 1.5) & (sampled[1000:] <= 2.0))


def test_cdf_rejects_empty_sampling_distribution():
    with pytest.raises(ValueError, match="no finite positive weight"):
        _build_cdf(np.array([0.0, 1.0]), np.array([-1.0, 1.0]), np.zeros((2, 2)))


def test_sample_charm_on_shell_and_in_bounds():
    rng = np.random.default_rng(1)
    pool = sample_meson_4vectors(5000, "charm", rng=rng)

    E = pool['E']
    px, py, pz = pool['px'], pool['py'], pool['pz']
    m_recon = np.sqrt(np.maximum(E**2 - (px**2 + py**2 + pz**2), 0.0))
    m_exp = np.array([MESON_MASSES[int(p)] for p in pool['species_pdg']])
    # Allow tiny numerical slop after bounded sub-bin smearing.
    assert np.allclose(m_recon, m_exp, atol=1e-5)

    # Kinematic envelope: smearing must not extend beyond the vendored table.
    assert pool['pt'].min() >= 0
    assert pool['pt'].max() <= 50.0
    assert abs(pool['y']).max() <= 3.0


def test_charm_species_distribution_matches_FRAG_C():
    rng = np.random.default_rng(2)
    pool = sample_meson_4vectors(20000, "charm", rng=rng)
    species = pool['species_pdg']
    expected = np.array([FRAG_C[421], FRAG_C[411], FRAG_C[431]])
    expected = expected / expected.sum()
    counts = np.array([(species == pdg).mean() for pdg in (421, 411, 431)])
    assert np.allclose(counts, expected, atol=0.02)


def test_bottom_species_distribution_matches_FRAG_B():
    rng = np.random.default_rng(3)
    pool = sample_meson_4vectors(20000, "bottom", rng=rng)
    species = pool['species_pdg']
    expected = np.array([FRAG_B[521], FRAG_B[511], FRAG_B[531]])
    expected = expected / expected.sum()
    counts = np.array([(species == pdg).mean() for pdg in (521, 511, 531)])
    assert np.allclose(counts, expected, atol=0.02)


def test_bc_pool_respects_bottom_table_bounds():
    pool = generate_bc_pool(5000, np.random.default_rng(4))
    assert pool['pt'].min() >= 0
    assert pool['pt'].max() <= 50.0
    assert abs(pool['y']).max() <= 3.0


def test_force_species_overrides_fragmentation_sampling():
    """force_species=541 must override species sampling and use M_BC for mass."""
    rng = np.random.default_rng(7)
    pool = sample_meson_4vectors(2000, "bottom", rng=rng, force_species=541)
    assert (pool['species_pdg'] == 541).all()
    E, px, py, pz = pool['E'], pool['px'], pool['py'], pool['pz']
    m_recon = np.sqrt(np.maximum(E**2 - (px**2 + py**2 + pz**2), 0.0))
    assert np.allclose(m_recon, MESON_MASSES[541], atol=1e-5)


def test_fragmentation_fractions_track_omitted_baryons():
    assert sum(FRAG_B.values()) < 1.0
    assert sum(FRAG_C.values()) < 1.0
    assert np.isclose(sum(FRAG_B.values()) + sum(OMITTED_FRAG_B.values()), 1.0)
    # ALICE Table 7 central values are published rounded; the meson + omitted
    # sum is 1.0007. atol=5e-3 absorbs that 7e-4 rounding gap with ~6x
    # headroom for digit-level updates to any single component. A genuine
    # ~20% relative shift in e.g. Xi_c+ (0.019 absolute) is deliberately
    # large enough to fail this test — that scale of change should not be
    # silently accommodated.
    assert np.isclose(
        sum(FRAG_C.values()) + sum(OMITTED_FRAG_C.values()), 1.0,
        rtol=0.0, atol=5e-3,
    )
