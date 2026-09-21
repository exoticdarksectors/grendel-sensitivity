"""Where the reference curves live, and the matplotlib environment."""
from __future__ import annotations

import os
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data"
HNL_DATA = DATA / "hnl"
BC4_DATA = DATA / "bc4"
BC10_DATA = DATA / "bc10"

# The HNLimits workbook (mhostert/Heavy-Neutrino-Limits) that carries the
# per-curve metadata the experimental bounds are transformed with.
HNLIMITS_DATABASE_URL = (
    "https://raw.githubusercontent.com/mhostert/Heavy-Neutrino-Limits/main/"
    "src/HNLimits/include/local_HNL_database.xlsx")
HNLIMITS_DATABASE = Path(os.environ.get(
    "GRENDEL_HNLIMITS_DATABASE", HNL_DATA / "local_HNL_database.xlsx")).expanduser()


def configure_matplotlib() -> None:
    """Select the non-interactive backend; called before pyplot is imported."""
    import matplotlib
    matplotlib.use("Agg")
