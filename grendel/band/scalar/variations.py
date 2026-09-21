"""The variations of the BC4 uncertainty campaign.

* one central FONLL grid;
* the six non-central members of the coherent seven-point scale set;
* 100 NNPDF4.0 NLO Monte-Carlo PDF replicas;
* two bottom-mass grids (4.5 and 5.0 GeV);
* one independently simulated decay-model alternate (LO-ChPT widths below
  2 GeV, perturbative spectator widths above) on the central grid;
* two same-physics central repeats with fresh seeds, the numerical controls,
  reported separately and never in the envelope.

The exHad decay model is a further ``decay_model``-axis member run by
``exhad_variation`` under the same independence policy.
"""
from __future__ import annotations

from pathlib import Path

from ..fonll_variations import assert_counts, bottom_records, select_variations  # noqa: F401

CENTRAL_SCHEME = "winkler"
DECAY_SCHEME = "chpt_spectator"
DECAY_VARIATION = "decay_chpt_spectator"
EXHAD_DECAY_SCHEME = "exhad:scalar-1809"
EXHAD_DECAY_VARIATION = "decay_exhad_1809"
NUMERICAL_CONTROL_VARIATIONS = ("central_repeat_1", "central_repeat_2")
EXPECTED_FONLL = {"central": 1, "scale": 6, "pdf": 100, "mass": 2}


def discover_variations(grid_dir: Path, validate_hashes=True) -> list[dict]:
    """The complete campaign, central first, the decay alternate second so
    the decay-model comparison completes before the long PDF campaign."""
    records = [{**r, "width_scheme": CENTRAL_SCHEME}
               for r in bottom_records(grid_dir, validate_hashes=validate_hashes)]
    assert_counts(records, EXPECTED_FONLL)
    central = records[0]
    decay = {**central, "name": DECAY_VARIATION, "axis": "decay_model", "width_scheme": DECAY_SCHEME}
    records = [central, decay, *records[1:]]
    records.extend({**central, "name": name, "axis": "numerical_control"}
                   for name in NUMERICAL_CONTROL_VARIATIONS)
    return records


def exhad_variation(grid_dir: Path, templates_dir: Path) -> dict:
    """The exHad decay-model member: the central grid with the exHad decay stage."""
    central = discover_variations(grid_dir)[0]
    return {**central, "name": EXHAD_DECAY_VARIATION, "axis": "decay_model",
            "width_scheme": EXHAD_DECAY_SCHEME, "decay_templates": Path(templates_dir)}
