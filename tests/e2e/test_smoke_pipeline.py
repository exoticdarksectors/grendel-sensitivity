"""The whole chain at tiny statistics for BC4, BC10 and the HNL."""
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


def _two_body_templates(path: Path, mass: float, ctau: float, pdg: tuple, masses: tuple, flavor: str) -> None:
    rng = np.random.default_rng(2)
    n_templates, (m1, m2) = 300, masses
    pstar = np.sqrt((mass**2 - (m1 + m2)**2) * (mass**2 - (m1 - m2)**2)) / (2 * mass)
    cos = rng.uniform(-1, 1, n_templates)
    azimuth = rng.uniform(0, 2 * np.pi, n_templates)
    sin = np.sqrt(1 - cos**2)
    d = np.column_stack([sin * np.cos(azimuth), sin * np.sin(azimuth), cos]) * pstar
    mom = np.empty((2 * n_templates, 3))
    mom[0::2], mom[1::2] = d, -d
    e = np.empty(2 * n_templates)
    e[0::2], e[1::2] = np.sqrt(pstar**2 + m1**2), np.sqrt(pstar**2 + m2**2)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, daughter_counts=np.full(n_templates, 2, np.int32),
                        pdg=np.tile(pdg, n_templates).astype(np.int32),
                        px=mom[:, 0], py=mom[:, 1], pz=mom[:, 2], energy=e,
                        mass=np.tile(masses, n_templates), charge=np.tile([-1.0, 1.0], n_templates),
                        stable=np.ones(2 * n_templates, bool), mass_GeV=np.float64(mass),
                        ctau_m_u2eq1=np.float64(ctau), n_templates=np.int32(n_templates), seed=np.int64(1),
                        flavor=np.array(flavor))


def _synthetic_hnl_inputs(root: Path, mass: float = 1.0) -> None:
    rng = np.random.default_rng(1)
    n = 4000
    p = rng.uniform(5.0, 60.0, n)
    eta = rng.uniform(-2.5, 2.5, n)
    phi = np.pi / 2 + rng.uniform(-0.6, 0.6, n)
    pt = p / np.cosh(eta)
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    energy = np.sqrt(p**2 + mass**2)
    weight = np.full(n, 1e8 / n)
    vectors = root / "vectors" / "Umu" / "combined"
    vectors.mkdir(parents=True)
    np.savetxt(vectors / "mN_1p000.csv", np.column_stack([weight, energy, px, py, pz]), delimiter=",", fmt="%.8e")
    _two_body_templates(root / "templates" / "Umu" / "templates_1p000.npz", mass, 1e-9, (13, 211),
                        (0.10566, 0.13957), "Umu")


def _synthetic_alp_templates(root: Path, masses) -> None:
    from grendel.io.vectors import format_mass_for_filename
    from grendel.models.alp import model
    for m in map(float, masses):
        _two_body_templates(root / f"templates_{format_mass_for_filename(m)}.npz", m,
                            model.alp_ctau(m, model.INV_F_REF), (13, -13), (0.10566, 0.10566), "BC10")


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
    _synthetic_alp_templates(work / "bc10" / "templates", MASSES)
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
    both = frame[frame["has_sensitivity_N10"]]
    assert len(both) > 0
    assert (both[f"{prefix}_min_N10"] >= both[f"{prefix}_min"] * (1 - 1e-12)).all()
    assert (both[f"{prefix}_max_N10"] <= both[f"{prefix}_max"] * (1 + 1e-12)).all()
