import numpy as np

from grendel.models.alp.production import (
    _importance_weighted_rate,
    _partition_parent_indices,
)


def test_alp_production_propagates_fonll_importance_weights():
    pool = {"sampling_weight": np.array([0.25, 0.5, 1.0, 2.0])}
    indices = np.array([3, 0, 2])

    weighted = _importance_weighted_rate(8.0, pool, indices)

    assert np.array_equal(weighted, np.array([16.0, 2.0, 8.0]))


def test_alp_production_defaults_to_nominal_unit_weights():
    weighted = _importance_weighted_rate(3.0, {}, np.array([2, 0, 1]))
    assert np.array_equal(weighted, np.full(3, 3.0))


def test_channel_parent_partition_has_no_bootstrap_duplicates():
    indices = _partition_parent_indices(
        101, 6, np.random.default_rng(12),
    )

    assert indices.shape == (6, 16)
    assert len(np.unique(indices)) == indices.size
    assert indices.min() >= 0
    assert indices.max() < 101
