"""Where the HNL production channels read and write.

The four-vector tree follows ``grendel.io.paths`` (``GRENDEL_HNL_VECTORS_DIR``
or ``<work>/hnl/llp_4vectors``); production adds a cache for the expensive
intermediate products:

    GRENDEL_HNL_CACHE_DIR         default <work>/hnl/cache
    GRENDEL_HNL_MG5_WORK_DIR      MadGraph run directories, default <cache>/madgraph
    GRENDEL_HNL_TAU_POOL_CSV      the prompt-tau pool, default <cache>/tau_pool.csv

Channel files are ``<vectors>/<flavor>/<channel>/mN_<label>.csv`` and the
per-point sum over channels is ``<vectors>/<flavor>/combined/mN_<label>.csv``.
"""
from __future__ import annotations

import os
from pathlib import Path

from ....io import vectors as _vectors
from ....io.paths import ModelPaths, work_dir
from ....io.vectors import (read_csv_matrix, scale_weight_column,  # noqa: F401
                            write_csv_matrix, write_empty_csv, write_llp_csv)


def _env_path(name: str, default: Path) -> Path:
    value = os.environ.get(name)
    return Path(value).expanduser() if value else default


LLP_VECTORS_DIR = ModelPaths.resolve("hnl").vectors
TMP_DIR = work_dir() / "hnl"
CACHE_DIR = _env_path("GRENDEL_HNL_CACHE_DIR", TMP_DIR / "cache")
MG5_WORK_DIR = _env_path("GRENDEL_HNL_MG5_WORK_DIR", CACHE_DIR / "madgraph")
TAU_POOL_CSV = _env_path("GRENDEL_HNL_TAU_POOL_CSV", CACHE_DIR / "tau_pool.csv")


def llp_csv_path(flavor, channel, mass, *, base=None, mkdir=True) -> Path:
    """``<base>/<flavor>/<channel>/mN_<label>.csv`` (base defaults to the
    four-vector tree)."""
    return _vectors.llp_csv_path(base or LLP_VECTORS_DIR, flavor, channel, mass, mkdir=mkdir)


def existing_tau_pool_csv() -> Path:
    """The prompt-tau pool, if one has been generated."""
    return TAU_POOL_CSV


def describe_paths() -> str:
    return (f"  four-vectors : {LLP_VECTORS_DIR}\n"
            f"  cache        : {CACHE_DIR}\n"
            f"  MadGraph work: {MG5_WORK_DIR}\n"
            f"  tau pool     : {TAU_POOL_CSV}")
