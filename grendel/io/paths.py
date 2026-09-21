"""Where a model reads its inputs and writes its outputs.

Event samples, four-vector pools, decay templates and ray-cast caches are far
too large for the repository. Every model resolves them through one policy:

    GRENDEL_WORK_DIR                     root for everything generated
                                         (default: ./grendel_work)
    GRENDEL_<MODEL>_VECTORS_DIR          four-vector CSVs, one per mass point
    GRENDEL_<MODEL>_TEMPLATES_DIR        rest-frame decay templates (.npz)
    GRENDEL_<MODEL>_ANALYSIS_DIR         scan output (sensitivity.csv, ...)
    GRENDEL_<MODEL>_GEOMETRY_DIR         ray-cast cache

with ``<MODEL>`` one of ``HNL``, ``BC4``, ``BC10``. Unset variables fall
back to ``<work>/<model>/{llp_4vectors,decay_templates,analysis,
analysis/geometry_cache}``. Command-line options override both.

The heavy-flavour production inputs have their own two:

    GRENDEL_FONLL_DIR                    a FONLL v1.3.3 source tree with the
                                         patched executables built
                                         (default: <third_party>/fonll)
    GRENDEL_FONLL_GRID_DIR               the FONLL variation campaign: grids
                                         plus variation_manifest.json
                                         (default: <work>/fonll_grids/output)

Paths are resolved when asked for, never at import time, so a process can
set the environment before it starts a scan.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

MODELS = ("hnl", "bc4", "bc10")


def work_dir() -> Path:
    return Path(os.environ.get("GRENDEL_WORK_DIR", "grendel_work")).expanduser()


def _env_path(name: str) -> Path | None:
    value = os.environ.get(name)
    return Path(value).expanduser() if value else None


@dataclass(frozen=True)
class ModelPaths:
    model: str
    vectors: Path
    templates: Path
    analysis: Path
    geometry: Path

    @classmethod
    def resolve(cls, model: str, *, vectors=None, templates=None,
                analysis=None, geometry=None) -> "ModelPaths":
        """Explicit arguments win, then the environment, then the work tree."""
        if model not in MODELS:
            raise ValueError(f"unknown model {model!r}; expected one of {MODELS}")
        key = model.upper()
        base = work_dir() / model
        analysis = Path(analysis) if analysis else \
            _env_path(f"GRENDEL_{key}_ANALYSIS_DIR") or base / "analysis"
        return cls(
            model=model,
            vectors=Path(vectors) if vectors else
            _env_path(f"GRENDEL_{key}_VECTORS_DIR") or base / "llp_4vectors",
            templates=Path(templates) if templates else
            _env_path(f"GRENDEL_{key}_TEMPLATES_DIR") or base / "decay_templates",
            analysis=analysis,
            geometry=Path(geometry) if geometry else
            _env_path(f"GRENDEL_{key}_GEOMETRY_DIR") or analysis / "geometry_cache",
        )

    def describe(self) -> str:
        return (f"{self.model}: vectors={self.vectors}\n"
                f"{' ' * len(self.model)}  templates={self.templates}\n"
                f"{' ' * len(self.model)}  analysis={self.analysis}\n"
                f"{' ' * len(self.model)}  geometry={self.geometry}")


def results_dir() -> Path | None:
    """Where finished scans are collected for the renderer
    (``GRENDEL_RESULTS_DIR``), if set."""
    return _env_path("GRENDEL_RESULTS_DIR")


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def third_party_dir() -> Path:
    """Where ``third_party/fetch.py`` checks out the pinned external packages
    (``GRENDEL_THIRD_PARTY_DIR`` or ``<repo>/third_party``)."""
    return _env_path("GRENDEL_THIRD_PARTY_DIR") or repo_root() / "third_party"


def fonll_dir() -> Path:
    """The FONLL source tree the grid generator drives (``GRENDEL_FONLL_DIR``
    or ``<third_party>/fonll``); see ``grendel.production.fonll_grids``."""
    return _env_path("GRENDEL_FONLL_DIR") or third_party_dir() / "fonll"


def fonll_grid_dir() -> Path:
    """Where the FONLL variation grids and their manifest live
    (``GRENDEL_FONLL_GRID_DIR`` or ``<work>/fonll_grids/output``)."""
    return _env_path("GRENDEL_FONLL_GRID_DIR") or work_dir() / "fonll_grids" / "output"
