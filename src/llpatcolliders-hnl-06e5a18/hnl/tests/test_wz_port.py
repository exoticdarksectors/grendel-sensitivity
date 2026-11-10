"""Data-path test for the W/Z -> l N port (no MadGraph required).

Builds a tiny synthetic unweighted LHE with HNL (PDG 9900012) final states,
runs it through LHEParser.write_hnl_csv, then applies the K-factor rescaling
exactly as run_wz_production.run_single_point does. Asserts:

  - the output schema is the shared 5-column weight,E,px,py,pz format,
  - the reconstructed invariant mass matches the injected HNL mass,
  - the summed weight scales by exactly K_FACTOR_EW.

This exercises the part of the port that has no external dependency, so it runs
in CI without the 148 MB MG5 install or Docker.
"""

import numpy as np
import pytest

from production.madgraph import runner
from production.madgraph.lhe_to_csv import LHEParser
from production.constants import K_FACTOR_EW

PDG_HNL = 9900012
M_N = 2.0
W_PER_EVENT = 3.5e-4  # XWGTUP = sigma_LO / N per unweighted event


def _hnl_line(px, py, pz):
    E = float(np.sqrt(px * px + py * py + pz * pz + M_N * M_N))
    # LHE particle line: id status mo1 mo2 col acol px py pz E m lifetime spin
    return (f" {PDG_HNL} 1 1 2 0 0 "
            f"{px:.9e} {py:.9e} {pz:.9e} {E:.9e} {M_N:.9e} 0. 1.")


def _write_synthetic_lhe(path, momenta):
    lines = ["<LesHouchesEvents version=\"3.0\">", "<init>", "</init>"]
    for (px, py, pz) in momenta:
        lines += [
            "<event>",
            # event header: nup idprup XWGTUP scalup aqedup aqcdup
            f" 1 1 {W_PER_EVENT:.9e} 9.118800e+01 7.546771e-03 1.300000e-01",
            _hnl_line(px, py, pz),
            "</event>",
        ]
    lines.append("</LesHouchesEvents>")
    path.write_text("\n".join(lines) + "\n")


def test_wz_lhe_to_csv_and_kfactor(tmp_path):
    momenta = [(1.0, 0.5, 10.0), (-2.0, 1.0, 30.0), (0.0, 0.0, 5.0)]
    lhe = tmp_path / "unweighted_events.lhe"
    _write_synthetic_lhe(lhe, momenta)

    csv_path = tmp_path / "mN_2p000.csv"
    n = LHEParser(lhe).write_hnl_csv(csv_path)
    assert n == len(momenta)

    data = np.loadtxt(csv_path, delimiter=",")
    assert data.ndim == 2 and data.shape == (len(momenta), 5)

    # Reconstructed invariant mass matches the injected HNL mass.
    E, px, py, pz = data[:, 1], data[:, 2], data[:, 3], data[:, 4]
    m_rec = np.sqrt(E**2 - px**2 - py**2 - pz**2)
    assert np.allclose(m_rec, M_N, atol=1e-6)

    # Pre-K-factor weights equal the LHE XWGTUP.
    assert np.allclose(data[:, 0], W_PER_EVENT, rtol=1e-9)
    w_sum_before = data[:, 0].sum()

    # Apply the K-factor exactly as run_single_point step 5 does.
    data[:, 0] *= K_FACTOR_EW
    np.savetxt(csv_path, data, delimiter=",", fmt="%.8e")

    rescaled = np.loadtxt(csv_path, delimiter=",")
    assert rescaled[:, 0].sum() == pytest.approx(w_sum_before * K_FACTOR_EW, rel=1e-6)


def test_wz_empty_lhe_yields_empty_csv(tmp_path):
    lhe = tmp_path / "unweighted_events.lhe"
    _write_synthetic_lhe(lhe, [])
    csv_path = tmp_path / "mN_2p000.csv"
    n = LHEParser(lhe).write_hnl_csv(csv_path)
    assert n == 0
    assert csv_path.exists() and csv_path.stat().st_size == 0


def test_mg5_process_generation_uses_cache_workdir(tmp_path, monkeypatch):
    work_dir = tmp_path / "cache"
    process_dir = work_dir / "process"
    proc_card = tmp_path / "proc_card.dat"
    proc_card.write_text("generate p p > w+\n")
    captured = {}

    def fake_run(*args, **kwargs):
        captured["cwd"] = kwargs["cwd"]
        (process_dir / "bin").mkdir(parents=True)
        (process_dir / "bin" / "generate_events").touch()
        return type("Result", (), {"returncode": 0})()

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    monkeypatch.setattr(runner, "force_compile_subprocesses", lambda *a, **k: True)
    monkeypatch.setattr(runner, "patch_rpath_for_lhapdf", lambda *a, **k: None)

    result = runner.ensure_process_dir(
        label="test",
        model_import="import model sm",
        proc_card=proc_card,
        work_dir=work_dir,
        process_dir=process_dir,
    )

    assert result == process_dir
    assert captured["cwd"] == work_dir
