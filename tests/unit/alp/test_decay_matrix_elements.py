import numpy as np

from grendel.models.alp.decay_matrix_elements import (
    CHANNEL_TO_MATRIX_ELEMENT,
    matrix_element_squared,
    normalized_template_weights,
)


def test_channel_mapping_covers_exported_three_body_decay_channels():
    assert CHANNEL_TO_MATRIX_ELEMENT["channel_005"] == "matrix_element_015"
    assert CHANNEL_TO_MATRIX_ELEMENT["channel_030"] == "matrix_element_011"
    assert len(CHANNEL_TO_MATRIX_ELEMENT) == 16


def test_cform_evaluation_matches_wolfram_anchor():
    value = matrix_element_squared("matrix_element_012", 2.0, [0.7], [0.5])
    assert np.isclose(value[0], 0.754344643799362, rtol=1e-13)


def test_weights_normalize_each_channel_without_changing_other_modes():
    channels = np.array([
        "channel_006", "channel_006", "channel_006", "channel_001"
    ])
    e1 = np.array([0.45, 0.55, 0.65, np.nan])
    e3 = np.array([0.35, 0.30, 0.25, np.nan])
    weights = normalized_template_weights(channels, e1, e3, 1.2)
    assert np.isclose(weights[:3].mean(), 1.0)
    assert weights[3] == 1.0
    assert np.all(weights >= 0.0)
