"""
Central path policy for HNL production.

Generated artifacts default to ``hnl/tmp`` so the source tree stays small:

  tmp/runs/<tag>/llp_4vectors/
  tmp/runs/<tag>/analysis/
  tmp/cache/madgraph/

Environment overrides:
  HNL_TMP_DIR, HNL_RUN_TAG, HNL_RUN_DIR, HNL_LLP_VECTORS_DIR,
  HNL_ANALYSIS_DIR, HNL_CACHE_DIR, HNL_MG5_WORK_DIR, HNL_TAU_POOL_CSV.
"""

from __future__ import annotations

import os
from pathlib import Path


HNL_ROOT = Path(__file__).resolve().parents[1]

RUN_TAG = os.environ.get("HNL_RUN_TAG", "default")

TMP_DIR = Path(os.environ.get("HNL_TMP_DIR", HNL_ROOT / "tmp")).expanduser()
RUN_DIR = Path(
    os.environ.get("HNL_RUN_DIR", TMP_DIR / "runs" / RUN_TAG)
).expanduser()
CACHE_DIR = Path(os.environ.get("HNL_CACHE_DIR", TMP_DIR / "cache")).expanduser()

LLP_VECTORS_DIR = Path(
    os.environ.get("HNL_LLP_VECTORS_DIR", RUN_DIR / "llp_4vectors")
).expanduser()
ANALYSIS_DIR = Path(
    os.environ.get("HNL_ANALYSIS_DIR", RUN_DIR / "analysis")
).expanduser()
MG5_WORK_DIR = Path(
    os.environ.get("HNL_MG5_WORK_DIR", CACHE_DIR / "madgraph")
).expanduser()
TAU_POOL_CSV = Path(
    os.environ.get("HNL_TAU_POOL_CSV", CACHE_DIR / "tau_pool.csv")
).expanduser()

# Compatibility with the previous committed/generated pool location. New
# generation writes to TAU_POOL_CSV; readers may fall back to this file.
LEGACY_TAU_POOL_CSV = HNL_ROOT / "vendored" / "tau_pool.csv"


def existing_tau_pool_csv() -> Path:
    """Return the preferred existing tau pool path, falling back to legacy."""
    if TAU_POOL_CSV.exists() and TAU_POOL_CSV.stat().st_size > 0:
        return TAU_POOL_CSV
    if LEGACY_TAU_POOL_CSV.exists() and LEGACY_TAU_POOL_CSV.stat().st_size > 0:
        return LEGACY_TAU_POOL_CSV
    return TAU_POOL_CSV


def describe_paths() -> str:
    """Human-readable summary for CLI logs."""
    return (
        f"run_dir={RUN_DIR}\n"
        f"llp_vectors={LLP_VECTORS_DIR}\n"
        f"analysis={ANALYSIS_DIR}\n"
        f"mg5_work={MG5_WORK_DIR}\n"
        f"tau_pool={TAU_POOL_CSV}"
    )
