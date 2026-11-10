"""End-to-end smoke: generate one mass point's CSV for one flavor/channel."""

import csv
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.external  # needs HNLCalc (third_party/fetch.py HNLCalc)


HNL_ROOT = Path(__file__).resolve().parent.parent




def test_generate_meson_csv_one_mass(tmp_path, monkeypatch):
    """Drive generate_meson_csvs.main() for Ue at m_N = 1.0 GeV."""
    monkeypatch.chdir(HNL_ROOT)
    monkeypatch.setenv("PYTHONUNBUFFERED", "1")

    # Redirect output to tmp_path so we don't clutter the repo.
    from grendel.models.hnl.production import generate_meson_csvs as G
    monkeypatch.setattr(G, "OUTPUT_BASE", tmp_path / "llp_4vectors")

    import sys
    argv_save = sys.argv
    try:
        sys.argv = [
            "generate_meson_csvs.py",
            "--flavor", "Ue",
            "--channel", "Dmeson",
            "--n-pool", "2000",
            "--masses", "1.0",
            "--seed", "0",
        ]
        G.main()
    finally:
        sys.argv = argv_save

    from grendel.io.vectors import format_mass_for_filename
    out = (tmp_path / "llp_4vectors" / "Ue" / "Dmeson"
           / f"mN_{format_mass_for_filename(1.0)}.csv")
    assert out.exists(), f"Expected CSV not produced at {out}"
    assert out.stat().st_size > 0, "CSV is empty; D -> N production should be open at 1.0 GeV"

    data = np.loadtxt(out, delimiter=",")
    assert data.ndim == 2 and data.shape[1] == 5
    weights = data[:, 0]
    assert (weights > 0).all()
    # Weighted sum should be a positive cross-section in pb, order-of-magnitude reasonable
    assert weights.sum() > 0
