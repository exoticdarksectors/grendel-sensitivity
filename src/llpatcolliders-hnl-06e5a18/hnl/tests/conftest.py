"""Shared pytest fixtures for hnl/tests."""

import sys
from pathlib import Path

# Make `hnl` importable as a package; tests address modules with full path.
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

# Make the hnl/ directory itself importable so existing modules can do
# `from config_mass_grid import ...` and `from production.X import ...`.
HNL_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HNL_ROOT))
sys.path.insert(0, str(HNL_ROOT / "vendored" / "HNLCalc"))
