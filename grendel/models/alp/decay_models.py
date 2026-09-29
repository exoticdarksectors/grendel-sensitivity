"""The pinned fermionic-ALP decay model behind BC10."""

from __future__ import annotations

from pathlib import Path


DEFAULT_DECAY_MODEL = "2501"
DECAY_MODEL_SPECS = {
    DEFAULT_DECAY_MODEL: {
        "data_subdir": "senscalc_2501",
        "phenomenology": "arXiv:2501.04525",
        "role": "central",
    },
}


def validate_decay_model(decay_model: str) -> str:
    if decay_model not in DECAY_MODEL_SPECS:
        choices = ", ".join(DECAY_MODEL_SPECS)
        raise ValueError(f"unknown decay model {decay_model!r}; choose {choices}")
    return decay_model


def decay_data_dir(decay_model: str = DEFAULT_DECAY_MODEL) -> Path:
    model = validate_decay_model(decay_model)
    return Path(__file__).resolve().parent / "data" / DECAY_MODEL_SPECS[model]["data_subdir"]
