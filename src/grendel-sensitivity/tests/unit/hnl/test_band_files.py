"""The HNL band tables are written where plot_money reads them: the FONLL
ribbon as hnl_band_fonll.csv, and the kaon d_esc band as kaon_desc_band.csv
with the three runs' lower edges in its columns."""
import json

import numpy as np
import pandas as pd
import pytest

from grendel.band import quadrature
from grendel.band.hnl import kaon_band
from grendel.models.hnl.plot_money import BAND_FILES


def _run(path, masses, u2_min, sensitive=None, open_=None):
    n = len(masses) * 2
    pd.DataFrame({
        "u2_min": np.tile(u2_min, 2), "u2_max": 1e-3, "has_sensitivity": True if sensitive is None
        else np.tile(sensitive, 2), "u2_min_open": False if open_ is None else np.tile(open_, 2),
        "u2_max_open": False, "mass_GeV": np.tile(masses, 2), "flavor": ["Umu"] * (n // 2) + ["Ue"] * (n // 2),
    }).to_csv(path, index=False)
    return path


def test_assembles_plot_money_columns(tmp_path):
    nominal = _run(tmp_path / "d15.csv", [0.2, 0.3, 1.0], [2e-8, 3e-8, 5e-9])
    near = _run(tmp_path / "d10.csv", [0.2, 0.3], [2.5e-8, 4e-8], sensitive=[True, False])
    far = _run(tmp_path / "d30.csv", [0.2, 0.3], [1.5e-8, 2e-8], open_=[False, True])
    out = tmp_path / "band" / BAND_FILES["kaon"]
    assert kaon_band.main(["--desc1p0", str(near), "--desc1p5", str(nominal), "--desc3p0", str(far),
                           "--out", str(out)]) == 0
    band = pd.read_csv(out)
    assert list(band.columns) == ["flavor", "mass_GeV", "u2_min_desc1p0", "u2_min_desc1p5", "u2_min_desc3p0"]
    assert list(band["flavor"]) == ["Ue", "Ue", "Umu", "Umu"]        # only the varied masses
    row = band[(band.flavor == "Umu") & (band.mass_GeV == 0.2)].iloc[0]
    assert (row.u2_min_desc1p0, row.u2_min_desc1p5, row.u2_min_desc3p0) == (2.5e-8, 2e-8, 1.5e-8)
    gap = band[(band.flavor == "Umu") & (band.mass_GeV == 0.3)].iloc[0]
    assert np.isnan(gap.u2_min_desc1p0) and np.isnan(gap.u2_min_desc3p0)   # no sensitivity / open edge
    assert gap.u2_min_desc1p5 == 3e-8


def test_mismatched_variation_points_fail(tmp_path):
    nominal = _run(tmp_path / "d15.csv", [0.2, 0.3], [2e-8, 3e-8])
    near = _run(tmp_path / "d10.csv", [0.2, 0.3], [2.5e-8, 4e-8])
    far = _run(tmp_path / "d30.csv", [0.2], [1.5e-8])
    with pytest.raises(ValueError, match="different"):
        kaon_band.assemble({1.0: near, 1.5: nominal, 3.0: far})
    short = _run(tmp_path / "d15s.csv", [0.2], [2e-8])
    with pytest.raises(ValueError, match="nominal"):
        kaon_band.assemble({1.0: near, 1.5: short, 3.0: near})


def test_quadrature_default_output_is_the_fonll_band_file(tmp_path):
    curve = _run(tmp_path / "central.csv", [1.0], [1e-8])
    registry = tmp_path / "band_registry.json"
    registry.write_text(json.dumps({"variations": [{"name": "central", "axis": "central",
                                                    "sensitivity_csv": str(curve)}]}))
    assert quadrature.main(["--registry", str(registry)]) == 0
    assert (tmp_path / BAND_FILES["fonll"]).is_file()
