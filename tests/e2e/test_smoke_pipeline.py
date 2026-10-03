"""The whole chain at tiny statistics for BC4, BC5, BC10 and the HNL."""
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

    _run(["grendel.models.scalar.scan_bc5", "--mass", "1.0", "20.0", "--decay-samples", "5",
          "--thresholds", "3", "10", "--n-pool", "3000", "--n-higgs", "3000",
          "--vectors-dir", str(work / "bc5" / "vectors"), "--out", str(out / "bc5")], env, repo)

    _run(["grendel.models.alp.production", "--mass", *MASSES, "--n-pool", "3000",
          "--out-dir", str(work / "bc10" / "vectors")], env, repo)
    _synthetic_alp_templates(work / "bc10" / "templates", MASSES)
    _run(["grendel.models.alp.scan", *common, "--vectors-dir", str(work / "bc10" / "vectors"),
          "--templates-dir", str(work / "bc10" / "templates"), "--out", str(out / "bc10")], env, repo)

    _synthetic_hnl_inputs(work / "hnl")
    _run(["grendel.models.hnl.scan", "--flavor", "Umu", "--mass", "1.0", "--decay-samples", "5",
          "--thresholds", "3", "10", "--vectors-dir", str(work / "hnl" / "vectors"),
          "--templates-dir", str(work / "hnl" / "templates"), "--out", str(out / "hnl")], env, repo)

    for model in ("hnl", "bc4", "bc5", "bc10"):
        _run(["grendel.io.thresholds", str(out / model / "sensitivity.csv"), "--threshold", "10"], env, repo)
    _run(["grendel.io.variants", str(out / "bc5" / "sensitivity.csv"), "--variant", "hSS"], env, repo)
    _run(["grendel.io.thresholds", str(out / "bc5" / "sensitivity_hSS.csv"), "--threshold", "10"], env, repo)
    return out, env, repo, work


@pytest.mark.parametrize("model,prefix", [("bc4", "u2"), ("bc5", "u2"), ("bc10", "invf"), ("hnl", "u2")])
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


BC5_VARIANTS = ("mixing", "hSS", "BSS", "brhss0.001", "brhss0.001_hSS", "brhss0.001_BSS")


def _same_edge(a, b, rtol=1e-5):
    return (np.isnan(a) and np.isnan(b)) or bool(np.isclose(a, b, rtol=rtol))


def test_bc5_splits_the_island_by_production_mode(results):
    out, _, _, _ = results
    frame = pd.read_csv(out / "bc5" / "sensitivity.csv").set_index("mass_GeV")
    for v in BC5_VARIANTS:
        assert {f"{v}_u2_min", f"{v}_u2_max", f"{v}_peak_N", f"{v}_peak_u2", f"{v}_has_sensitivity",
                f"{v}_u2_min_open", f"{v}_u2_max_open", f"{v}_u2_min_N10", f"{v}_has_sensitivity_N10"} <= set(frame.columns)
    assert {"n_hits_mixing", "n_hits_hSS", "n_hits_BSS", "n_hits_quartic",
            "u2_min_N_mixing", "u2_min_N_quartic", "u2_min_N_hSS", "peak_u2_N_BSS"} <= set(frame.columns)
    for label in ("u2_min", "peak_u2", "u2_max"):
        quartic, modes = frame[f"{label}_N_quartic"], frame[f"{label}_N_hSS"] + frame[f"{label}_N_BSS"]
        known = quartic.notna()
        assert known.any() and np.allclose(quartic[known], modes[known], rtol=1e-9)
    assert (frame["n_hits_mixing"] + frame["n_hits_hSS"] + frame["n_hits_BSS"] == frame["n_hits_mixing"] + frame["n_hits_quartic"]).all()

    # every curve solved on a subset (or a down-scaled copy) of the nominal rows lies inside the nominal island
    for mass, row in frame.iterrows():
        for v in BC5_VARIANTS:
            assert row[f"{v}_peak_N"] <= row["peak_N"] * (1 + 1e-6), (mass, v)
            if not row[f"{v}_has_sensitivity"]:
                continue
            if np.isfinite(row["u2_min"]):
                assert row[f"{v}_u2_min"] >= row["u2_min"] * (1 - 1e-5), (mass, v)
            if np.isfinite(row["u2_max"]):
                assert row[f"{v}_u2_max"] <= row["u2_max"] * (1 + 1e-5), (mass, v)
    assert np.allclose(frame["brhss0.001_hSS_peak_N"], 0.1 * frame["hSS_peak_N"], rtol=1e-6)

    # above m_B only the Higgs makes S: no mixing or B -> X SS rows, and the h -> SS island is the island
    high = frame.loc[20.0]
    assert high["n_hits_mixing"] == 0 and high["n_hits_BSS"] == 0 and high["n_hits_hSS"] > 0
    assert not high["mixing_has_sensitivity"] and not high["BSS_has_sensitivity"]
    assert high["hSS_has_sensitivity"] and high["has_sensitivity"]
    for field in ("u2_min", "u2_max", "peak_u2", "peak_N"):
        assert _same_edge(high[f"hSS_{field}"], high[field]), field
    # at 1 GeV every mode produces S
    low = frame.loc[1.0]
    assert low["n_hits_mixing"] > 0 and low["n_hits_hSS"] > 0 and low["n_hits_BSS"] > 0


def test_bc5_variant_projection_has_the_nominal_layout(results):
    out, _, _, _ = results
    root = out / "bc5"
    frame = pd.read_csv(root / "sensitivity.csv")
    hss = pd.read_csv(root / "sensitivity_hSS.csv")
    assert list(hss["mass_GeV"]) == list(frame["mass_GeV"])
    assert {"u2_min", "u2_max", "peak_N", "peak_u2", "has_sensitivity", "u2_min_open", "u2_min_N10",
            "has_sensitivity_N10", "n_hits", "u2_min_N_mixing"} <= set(hss.columns)
    assert not any(c.startswith(("hSS_", "mixing_", "BSS_", "brhss")) for c in hss.columns)
    for field in ("u2_min", "u2_max", "peak_N", "u2_min_N10"):
        assert np.allclose(hss[field], frame[f"hSS_{field}"], equal_nan=True), field
    nsig10 = pd.read_csv(root / "sensitivity_hSS_nsig10.csv")
    assert "u2_min_N10" not in nsig10.columns
    assert np.allclose(nsig10["u2_min"], frame["hSS_u2_min_N10"], equal_nan=True)
