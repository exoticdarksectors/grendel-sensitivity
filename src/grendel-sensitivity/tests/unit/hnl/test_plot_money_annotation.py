"""The N -> mu pi closing step (m_mu + m_pi = 0.245 GeV) is annotated on the
Umu panel only: tau mixing has no such step."""
import matplotlib.figure
import pandas as pd

from grendel.models.hnl import plot_money


def test_pion_step_annotation_is_on_the_umu_panel_only(tmp_path, monkeypatch):
    central = tmp_path / "hnl" / "sensitivity.csv"
    central.parent.mkdir()
    pd.DataFrame([{"flavor": flavor, "mass_GeV": mass, "has_sensitivity": True, "u2_min": 1e-8,
                   "u2_max": 1e-4, "u2_min_open": False, "u2_max_open": False}
                  for flavor in ("Ue", "Umu", "Utau") for mass in (0.2, 0.3, 1.0, 2.0)]).to_csv(central, index=False)
    figures = []
    monkeypatch.setattr(matplotlib.figure.Figure, "savefig", lambda self, *args, **kwargs: figures.append(self))
    assert plot_money.main(["--central-csv", str(central), "--out-dir", str(tmp_path / "out")]) == 0
    annotated = [ax.get_title() for ax in figures[0].axes if any("closes" in t.get_text() for t in ax.texts)]
    assert annotated == [f"HNL {plot_money.LAB['Umu']}"]
