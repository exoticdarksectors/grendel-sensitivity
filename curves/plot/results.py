"""The scan outputs a figure is drawn from.

``python -m grendel.models.<model>.scan --out RESULTS/<model>`` leaves
``sensitivity.csv`` there; the same scan solved at a second threshold is
``sensitivity_nsig10.csv``, and the band campaigns write under ``band/``::

    RESULTS/
      hnl/sensitivity.csv    hnl/sensitivity_nsig10.csv    hnl/band/...
      bc4/sensitivity.csv    bc4/sensitivity_nsig10.csv    bc4/band/...
      bc10/sensitivity.csv   bc10/sensitivity_nsig10.csv   bc10/band/...

The renderer takes the root as ``--grendel-dir`` or ``GRENDEL_RESULTS_DIR``.
"""
from __future__ import annotations

import os
from pathlib import Path

MODELS = ("hnl", "bc4", "bc10")
ENVELOPE_FILES = {
    "bc4": "bc4_single_source_variation_envelope.csv",
    "bc10": "bc10_single_source_variation_envelope.csv",
}
# The HNL variation diagnostics (grendel.band.hnl): the FONLL production
# ribbon, the decay-model legs and the B_c normalisation nuisance.
HNL_BAND_FILES = {
    "fonll": "hnl_band_fonll.csv",
    "decay": "decay_model_band.csv",
    "bc": "bc_nuisance_band.csv",
}


class GrendelResults:
    def __init__(self, root=None):
        root = root or os.environ.get("GRENDEL_RESULTS_DIR")
        if not root:
            raise FileNotFoundError(
                "no results directory: pass --grendel-dir or set GRENDEL_RESULTS_DIR "
                "to the tree holding <model>/sensitivity.csv")
        self.root = Path(root).expanduser()

    def _require(self, path: Path, what: str) -> Path:
        if not path.is_file():
            raise FileNotFoundError(f"{what} not found: {path}")
        return path

    def sensitivity(self, model: str) -> Path:
        return self._require(self.root / model / "sensitivity.csv", f"{model} sensitivity curve")

    def nsig10(self, model: str) -> Path | None:
        """The same scan solved at N_sig >= 10, if present."""
        path = self.root / model / "sensitivity_nsig10.csv"
        return path if path.is_file() else None

    def envelope(self, model: str) -> Path:
        return self._require(self.root / model / "band" / ENVELOPE_FILES[model],
                             f"{model} variation envelope")

    def band(self, model: str, name: str) -> Path:
        return self._require(self.root / model / "band" / name, f"{model} band {name}")

    def describe(self) -> str:
        lines = [f"results: {self.root}"]
        for model in MODELS:
            s = self.root / model / "sensitivity.csv"
            lines.append(f"  {model:4s} {'ok' if s.is_file() else 'missing'}  {s}")
        return "\n".join(lines)
