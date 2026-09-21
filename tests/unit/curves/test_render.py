"""Rendering from a synthetic results tree."""
import csv
from pathlib import Path

import pytest

from curves.plot.results import GrendelResults

BC_FIELDS = ["mass_GeV", "has_sensitivity", "peak_N", "{p}_min", "{p}_max", "{p}_min_open", "{p}_max_open",
             "peak_{p}"]


def _write_bc(path: Path, prefix: str, scale: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [f.format(p=prefix) for f in BC_FIELDS]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for i, mass in enumerate((0.5, 0.75, 1.0, 1.5, 2.0, 2.5)):
            sensitive = i < 5
            peak = 30.0 / (1 + i) if sensitive else 2.0
            w.writerow({"mass_GeV": mass, "has_sensitivity": sensitive, "peak_N": peak,
                        f"{prefix}_min": 1e-9 * scale * (1 + i) if sensitive else "",
                        f"{prefix}_max": 1e-6 / scale / (1 + i) if sensitive else "",
                        f"{prefix}_min_open": False, f"{prefix}_max_open": False,
                        f"peak_{prefix}": 1e-8 * (1 + i)})


@pytest.fixture
def results(tmp_path):
    for model, prefix in (("bc4", "u2"), ("bc10", "invf")):
        _write_bc(tmp_path / model / "sensitivity.csv", prefix, 1.0)
        _write_bc(tmp_path / model / "sensitivity_nsig10.csv", prefix, 1.8)
    return GrendelResults(tmp_path)


def test_results_layout_is_explicit(tmp_path):
    with pytest.raises(FileNotFoundError, match="sensitivity curve"):
        GrendelResults(tmp_path).sensitivity("bc4")
    assert GrendelResults(tmp_path).nsig10("bc4") is None


@pytest.mark.parametrize("benchmark", ["bc4", "bc10"])
def test_paper_figure_renders_from_a_results_tree(results, tmp_path, benchmark):
    from curves.plot.paper import render_bc_paper
    out = render_bc_paper(benchmark, results, tmp_path / "figures")
    assert out.with_suffix(".png").stat().st_size > 0
    assert out.stat().st_size > 0


def test_talk_figure_renders_from_a_results_tree(results, tmp_path):
    from curves.plot.talk import render_bc_talk
    out = render_bc_talk("bc4", results, tmp_path / "figures")
    assert out.with_suffix(".png").stat().st_size > 0


def test_cli_paper_layout(results, tmp_path):
    from curves.plot.__main__ import main
    assert main(["paper", "--grendel-dir", str(results.root), "--out", str(tmp_path / "f"),
                 "--benchmark", "bc4", "--suffix", "_x"]) == 0
    assert (tmp_path / "f" / "bc4_grendel_paper_x.png").is_file()


def test_hnl_panel_renders_with_hnlimits(results, tmp_path):
    pytest.importorskip("HNLimits")
    from curves.plot.paper import render_hnl_paper
    path = results.root / "hnl" / "sensitivity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["flavor", "mass_GeV", "u2_min", "u2_max", "u2_min_open", "u2_max_open", "peak_N", "peak_u2",
              "has_sensitivity"]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for flavor in ("Ue", "Umu", "Utau"):
            for i, mass in enumerate((0.5, 1.0, 2.0, 3.0)):
                sensitive = i < 3
                w.writerow({"flavor": flavor, "mass_GeV": mass, "u2_min": 1e-8 * (1 + i) if sensitive else "",
                            "u2_max": 1e-5 if sensitive else "", "u2_min_open": False, "u2_max_open": False,
                            "peak_N": 20.0 / (1 + i) if sensitive else 1.0, "peak_u2": 1e-7,
                            "has_sensitivity": sensitive})
    out = render_hnl_paper(results, tmp_path / "figures")
    assert out.with_suffix(".png").stat().st_size > 0
