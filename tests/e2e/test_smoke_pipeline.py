"""The whole chain at tiny statistics: produce -> scan -> split thresholds
-> render, for all three models, through the command-line entry points.

BC4 and BC10 produce their four-vectors from the tracked FONLL grid; BC10
uses the two-track proxy templates. The HNL production needs HNLCalc, so
its four-vectors and rest-frame templates are synthesised here: a beam of
heavy neutral leptons aimed at the detector decaying to mu pi. The numbers
mean nothing physically; what is checked is that every stage runs, writes
the documented layout, and that the N >= 10 island sits inside the
N >= 3 one.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

MASSES = ["1.0", "2.0"]


def _run(args, env, cwd):
    result = subprocess.run([sys.executable, "-m", *args], env=env, cwd=cwd, capture_output=True, text=True)
    assert result.returncode == 0, f"{' '.join(args)} failed:\n{result.stdout[-2000:]}\n{result.stderr[-3000:]}"
    return result


def _synthetic_hnl_inputs(root: Path, mass: float = 1.0) -> None:
    rng = np.random.default_rng(1)
    n = 4000
    p = rng.uniform(5.0, 60.0, n)
    eta = rng.uniform(-2.5, 2.5, n)
    phi = np.pi / 2 + rng.uniform(-0.6, 0.6, n)          # towards the gallery above the IP
    pt = p / np.cosh(eta)
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    energy = np.sqrt(p**2 + mass**2)
    weight = np.full(n, 1e8 / n)                          # pb per event
    vectors = root / "vectors" / "Umu" / "combined"
    vectors.mkdir(parents=True)
    np.savetxt(vectors / "mN_1p000.csv", np.column_stack([weight, energy, px, py, pz]), delimiter=",", fmt="%.8e")

    n_templates, m_mu, m_pi = 300, 0.10566, 0.13957
    pstar = np.sqrt((mass**2 - (m_mu + m_pi)**2) * (mass**2 - (m_mu - m_pi)**2)) / (2 * mass)
    cos = rng.uniform(-1, 1, n_templates)
    azimuth = rng.uniform(0, 2 * np.pi, n_templates)
    sin = np.sqrt(1 - cos**2)
    d = np.column_stack([sin * np.cos(azimuth), sin * np.sin(azimuth), cos]) * pstar
    # daughters interleaved per decay: mu then pi
    mom = np.empty((2 * n_templates, 3))
    mom[0::2], mom[1::2] = d, -d
    e = np.empty(2 * n_templates)
    e[0::2], e[1::2] = np.sqrt(pstar**2 + m_mu**2), np.sqrt(pstar**2 + m_pi**2)
    templates = root / "templates" / "Umu"
    templates.mkdir(parents=True)
    np.savez_compressed(templates / "templates_1p000.npz",
                        daughter_counts=np.full(n_templates, 2, np.int32),
                        pdg=np.tile([13, 211], n_templates).astype(np.int32),
                        px=mom[:, 0], py=mom[:, 1], pz=mom[:, 2], energy=e,
                        mass=np.tile([m_mu, m_pi], n_templates), charge=np.tile([-1.0, 1.0], n_templates),
                        stable=np.ones(2 * n_templates, bool), mass_GeV=np.float64(mass),
                        ctau_m_u2eq1=np.float64(1e-9), n_templates=np.int32(n_templates), seed=np.int64(1),
                        flavor=np.array("Umu"))


@pytest.fixture(scope="module")
def results(tmp_path_factory):
    work = tmp_path_factory.mktemp("smoke")
    env = {**os.environ, "GRENDEL_WORK_DIR": str(work / "work"), "OMP_NUM_THREADS": "1",
           "MPLCONFIGDIR": str(work / "mpl")}
    repo = Path(__file__).resolve().parents[2]
    out = work / "results"
    common = ["--mass", *MASSES, "--decay-samples", "5", "--thresholds", "3", "10"]

    _run(["grendel.models.scalar.scan", *common, "--n-pool", "3000",
          "--vectors-dir", str(work / "bc4" / "vectors"), "--out", str(out / "bc4")], env, repo)

    _run(["grendel.models.alp.production", "--mass", *MASSES, "--n-pool", "3000",
          "--out-dir", str(work / "bc10" / "vectors")], env, repo)
    _run(["grendel.models.alp.templates_proxy", "--legacy-proxy", "--mass", *MASSES, "--n-templates", "300",
          "--out", str(work / "bc10" / "templates")], env, repo)
    _run(["grendel.models.alp.scan", *common, "--vectors-dir", str(work / "bc10" / "vectors"),
          "--templates-dir", str(work / "bc10" / "templates"), "--out", str(out / "bc10")], env, repo)

    _synthetic_hnl_inputs(work / "hnl")
    _run(["grendel.models.hnl.scan", "--flavor", "Umu", "--mass", "1.0", "--decay-samples", "5",
          "--thresholds", "3", "10", "--vectors-dir", str(work / "hnl" / "vectors"),
          "--templates-dir", str(work / "hnl" / "templates"), "--out", str(out / "hnl")], env, repo)

    for model in ("hnl", "bc4", "bc10"):
        _run(["grendel.io.thresholds", str(out / model / "sensitivity.csv"), "--threshold", "10"], env, repo)
    return out, env, repo, work


@pytest.mark.parametrize("model,prefix", [("bc4", "u2"), ("bc10", "invf"), ("hnl", "u2")])
def test_scan_writes_the_documented_layout(results, model, prefix):
    out, _, _, _ = results
    root = out / model
    assert {p.name for p in root.iterdir()} >= {"sensitivity.csv", "sensitivity_nsig10.csv", "run_metadata.json",
                                                "scan_status.json", "geometry_cache"}
    assert list((root / "geometry_cache").rglob("*.npz"))
    meta = json.loads((root / "run_metadata.json").read_text())
    assert meta["model"] == model and meta["thresholds"] == [3.0, 10.0]
    frame = pd.read_csv(root / "sensitivity.csv")
    assert {"mass_GeV", "has_sensitivity", "peak_N", f"{prefix}_min", f"{prefix}_max", f"{prefix}_min_open",
            f"{prefix}_max_open", f"{prefix}_min_N10", "has_sensitivity_N10"} <= set(frame.columns)
    assert frame["has_sensitivity"].all(), frame[["mass_GeV", "peak_N"]]
    nsig10 = pd.read_csv(root / "sensitivity_nsig10.csv")
    assert f"{prefix}_min_N10" not in nsig10.columns and list(nsig10["mass_GeV"]) == list(frame["mass_GeV"])
    # The stricter island is nested inside the baseline one.
    both = frame[frame["has_sensitivity_N10"]]
    assert len(both) > 0
    assert (both[f"{prefix}_min_N10"] >= both[f"{prefix}_min"] * (1 - 1e-12)).all()
    assert (both[f"{prefix}_max_N10"] <= both[f"{prefix}_max"] * (1 + 1e-12)).all()


def test_paper_figures_render_from_the_results_tree(results, tmp_path):
    out, env, repo, _ = results
    benchmarks = ["bc4", "bc10"]
    try:
        import HNLimits  # noqa: F401
        benchmarks.append("hnl")
    except ImportError:
        pass
    _run(["curves.plot", "paper", "--grendel-dir", str(out), "--out", str(tmp_path), "--benchmark", *benchmarks],
         env, repo)
    for benchmark in benchmarks:
        stem = "hnlimits_grendel_paper" if benchmark == "hnl" else f"{benchmark}_grendel_paper"
        assert (tmp_path / f"{stem}.png").stat().st_size > 0
        assert (tmp_path / f"{stem}.pdf").stat().st_size > 0


def test_cutflow_reproduces_the_scan_selection(results, tmp_path):
    out, env, repo, work = results
    _run(["grendel.cutflow", "--model", "bc4", "--mass", "1.0", "--from-results", str(out / "bc4" / "sensitivity.csv"),
          "--decay-samples", "5", "--vectors-dir", str(work / "bc4" / "vectors"), "--out", str(out / "bc4"),
          "--table", str(tmp_path / "cutflow.csv")], env, repo)
    table = pd.read_csv(tmp_path / "cutflow.csv")
    assert len(table) > 0
