"""Evaluate exact three-body matrix elements from pinned SensCalc models."""

from __future__ import annotations

import json
import sys
from functools import lru_cache

import numpy as np

from ...io.tables import read_table_text
from .decay_models import DEFAULT_DECAY_MODEL, decay_data_dir


DATA_PATH = (
    decay_data_dir(DEFAULT_DECAY_MODEL) / "matrix_elements.json"
)

CHANNEL_TO_MATRIX_ELEMENT = {
    "channel_005": "matrix_element_015",
    "channel_006": "matrix_element_012",
    "channel_009": "matrix_element_013",
    "channel_010": "matrix_element_014",
    "channel_013": "matrix_element_016",
    "channel_014": "matrix_element_004",
    "channel_015": "matrix_element_017",
    "channel_016": "matrix_element_018",
    "channel_021": "matrix_element_001",
    "channel_022": "matrix_element_002",
    "channel_023": "matrix_element_006",
    "channel_025": "matrix_element_007",
    "channel_026": "matrix_element_008",
    "channel_027": "matrix_element_009",
    "channel_028": "matrix_element_010",
    "channel_030": "matrix_element_011",
}


def _unit_step(value):
    return np.heaviside(value, 1.0)


_EVAL_GLOBALS = {
    "__builtins__": {},
    "Power": np.power,
    "Sqrt": np.sqrt,
    "Complex": complex,
    "UnitStep": _unit_step,
}


@lru_cache(maxsize=None)
def _compiled_expressions(
    decay_model: str = DEFAULT_DECAY_MODEL,
) -> dict[str, object]:
    records = json.loads(
        read_table_text(decay_data_dir(decay_model) / "matrix_elements.json")
    )
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 200_000))
    return {
        record["id"]: compile(
            record["expression_c_form"], record["id"], "eval"
        )
        for record in records
    }


def matrix_element_squared(
    matrix_element_id: str,
    mass_gev: float,
    energy_1_gev,
    energy_3_gev,
    decay_model: str = DEFAULT_DECAY_MODEL,
) -> np.ndarray:
    """Return the real, non-negative squared amplitude on a Dalitz sample."""
    values = eval(
        _compiled_expressions(decay_model)[matrix_element_id],
        _EVAL_GLOBALS,
        {
            "mLLP": float(mass_gev),
            "E1": np.asarray(energy_1_gev, dtype=float),
            "E3": np.asarray(energy_3_gev, dtype=float),
        },
    )
    values = np.real(np.asarray(values, dtype=complex)).astype(float)
    return np.where(np.isfinite(values), np.maximum(values, 0.0), 0.0)


def normalized_template_weights(
    channel_ids,
    energy_1_gev,
    energy_3_gev,
    mass_gev: float,
    decay_model: str = DEFAULT_DECAY_MODEL,
) -> np.ndarray:
    """Matrix-element reweights normalized within every exclusive channel."""
    channel_ids = np.asarray(channel_ids).astype(str)
    energy_1_gev = np.asarray(energy_1_gev, dtype=float)
    energy_3_gev = np.asarray(energy_3_gev, dtype=float)
    if not (channel_ids.shape == energy_1_gev.shape == energy_3_gev.shape):
        raise ValueError("channel and primary-energy arrays must have equal shape")

    weights = np.ones(channel_ids.shape, dtype=float)
    for channel_id, matrix_id in CHANNEL_TO_MATRIX_ELEMENT.items():
        selected = channel_ids == channel_id
        if not selected.any():
            continue
        raw = matrix_element_squared(
            matrix_id,
            mass_gev,
            energy_1_gev[selected],
            energy_3_gev[selected],
            decay_model,
        )
        mean = float(raw.mean())
        if mean > 0.0:
            weights[selected] = raw / mean
    return weights
