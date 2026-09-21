"""Machine-readable exclusive definitions for the pinned ALP decay models."""

from __future__ import annotations

import json
from functools import lru_cache

import numpy as np

from ...io.tables import open_table
from .decay_models import DEFAULT_DECAY_MODEL, decay_data_dir


DATA_DIR = decay_data_dir(DEFAULT_DECAY_MODEL)
HBAR_C_GEV_M = 1.973269804e-16
INV_F_REF = 1.0e-3

# Stable or Pythia-known primary PDG products corresponding positionally to
# decay_channels.json. Pythia subsequently decays unstable mesons and taus and
# showers/hadronizes quark pairs. The isolated two-gluon channel is replaced
# by the configured light-quark surrogate during template generation.
DECAY_PRODUCTS_PDG = {
    "channel_001": (-11, 11),
    "channel_002": (-13, 13),
    "channel_003": (-15, 15),
    "channel_004": (22, 22),
    "channel_005": (211, -211, 111),
    "channel_006": (211, -211, 22),
    "channel_007": (211, 111, -211, 111),
    "channel_008": (211, -211, 211, -211),
    "channel_009": (221, 111, 111),
    "channel_010": (221, 211, -211),
    "channel_011": (323, -323),
    "channel_012": (313, -313),
    "channel_013": (223, 211, -211),
    "channel_014": (111, 111, 111),
    "channel_015": (331, 111, 111),
    "channel_016": (331, 211, -211),
    "channel_017": (223, 223),
    "channel_018": (21, 21),
    "channel_019": (4, -4),
    "channel_020": (3, -3),
    "channel_021": (130, 130, 111),
    "channel_022": (310, 310, 111),
    "channel_023": (130, 310, 111),
    "channel_024": (313, -313),
    "channel_025": (-321, 130, 211),
    "channel_026": (321, 130, -211),
    "channel_027": (-321, 310, 211),
    "channel_028": (321, 310, -211),
    "channel_029": (323, -323),
    "channel_030": (321, -321, 111),
    "channel_031": (113, 113),
    "channel_032": (213, -213),
}

DUPLICATE_CHANNEL_IDS = {"channel_024", "channel_029"}
UNCLASSIFIED_CHANNEL_ID = "unclassified_neutral_remainder"

@lru_cache(maxsize=None)
def _branching_table(decay_model: str = DEFAULT_DECAY_MODEL):
    return np.genfromtxt(
        open_table(decay_data_dir(decay_model) / "branching_ratios.csv"),
        delimiter=",",
        names=True,
    )


@lru_cache(maxsize=None)
def _total_width_table(decay_model: str = DEFAULT_DECAY_MODEL):
    data_dir = decay_data_dir(decay_model)
    metadata = json.loads((data_dir / "widths_metadata.json").read_text())
    total = next(
        entry for entry in metadata["columns"]
        if entry["canonical_name"] == "total"
    )
    table = np.loadtxt(open_table(data_dir / "widths_bnt.csv"), delimiter=",", skiprows=1)
    return table[:, 0], table[:, total["source_index"] - 1]


def exclusive_branching_weights(
    mass_gev: float, decay_model: str = DEFAULT_DECAY_MODEL
) -> dict[str, float]:
    """Unique exclusive BRs plus a conservative neutral missing-width mode."""
    table = _branching_table(decay_model)
    masses = table["mass_GeV"]
    if not masses[0] <= mass_gev <= masses[-1]:
        raise ValueError(f"mass {mass_gev:g} GeV is outside the decay table")
    weights = {
        channel_id: max(float(np.interp(mass_gev, masses, table[channel_id])), 0.0)
        for channel_id in DECAY_PRODUCTS_PDG
        if channel_id not in DUPLICATE_CHANNEL_IDS
    }
    known = sum(weights.values())
    if known > 1.0:
        weights = {key: value / known for key, value in weights.items()}
        known = 1.0
    weights[UNCLASSIFIED_CHANNEL_ID] = max(1.0 - known, 0.0)
    return {key: value for key, value in weights.items() if value > 0.0}


def ctau_at_reference_coupling(
    mass_gev: float, decay_model: str = DEFAULT_DECAY_MODEL
) -> float:
    masses, coefficients = _total_width_table(decay_model)
    if not masses[0] <= mass_gev <= masses[-1]:
        raise ValueError(f"mass {mass_gev:g} GeV is outside the width table")
    width = float(np.interp(mass_gev, masses, coefficients)) * INV_F_REF**2
    return HBAR_C_GEV_M / width if width > 0.0 else np.inf
