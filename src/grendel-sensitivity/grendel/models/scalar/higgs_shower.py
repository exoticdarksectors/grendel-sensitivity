"""BC5 Higgs production, stage 2: shower the MadGraph LHE events with
Pythia 8 and record the Higgs four-vector after the shower.

    python -m grendel.models.scalar.higgs_shower            # every mode in the manifest
    python -m grendel.models.scalar.higgs_shower --modes ggf --max-events 2000

Reads ``higgs_lhe_manifest.json`` of ``higgs_production`` (in
``$GRENDEL_BC5_MG5_WORK_DIR``), runs the parton shower on each LHE with the
Higgs kept stable (ISR + FSR, which give the Higgs its recoil; no
hadronisation and no MPI, which do not move it), and writes ``higgs_<mode>.csv`` (columns
``pt_GeV, y, E, px, py, pz``) next to the manifest, plus ``higgs_shower_manifest.json``
with the Pythia version, the event counts and the per-mode cross section
Pythia reports from the LHE. The shower is ``higgs_shower.cc``, compiled
against Pythia 8.317 (the version exHad needs; ``pythia8-config`` from
``$GRENDEL_PYTHIA8_CONFIG`` or ``<third_party>/pythia8317``, as for the kaon
spectrum), so no Pythia Python bindings are needed.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path

import numpy as np

from ..hnl.production.make_kaon_spectrum import _cfg, pythia8_config
from .higgs_production import mg5_work_dir

HIGGS_PDG = 25
SRC = Path(__file__).with_suffix(".cc")


def _plain_lhe(lhe: Path, tmp_dir: Path) -> Path:
    """Pythia's LHEF reader wants a plain file unless built with gzip."""
    if lhe.suffix != ".gz":
        return lhe
    out = tmp_dir / lhe.with_suffix("").name
    with gzip.open(lhe, "rb") as src, open(out, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return out


def shower_lhe(lhe: Path, *, seed: int = 1, max_events: int | None = None, quiet: bool = True):
    cfg = pythia8_config()
    with tempfile.TemporaryDirectory() as tmp:
        plain = _plain_lhe(lhe, Path(tmp))
        binary = Path(tmp) / "higgs_shower"
        subprocess.run(["c++", "-O2", str(SRC), "-o", str(binary), *_cfg(cfg, "--cxxflags"), *_cfg(cfg, "--libs")],
                       check=True)
        env = dict(os.environ, PYTHIA8DATA=subprocess.check_output([cfg, "--xmldoc"], text=True).strip())
        proc = subprocess.run([str(binary), str(plain), str(seed), str(max_events or 0)],
                              capture_output=True, text=True, env=env, check=True)
    meta = dict(line.split() for line in proc.stderr.splitlines() if len(line.split()) == 2 and line.split()[0].isupper())
    rows = np.array([line.split()[1:] for line in proc.stdout.splitlines() if line.startswith("H ")], float)
    summary = {
        "pythia_version": float(meta["VERSION"]),
        "n_showered": int(meta["NSHOWERED"]),
        "n_higgs": int(meta["NHIGGS"]),
        "sigma_lhe_pb": float(meta["SIGMAGEN_MB"]) * 1e9,     # Pythia reports mb
        "seed": seed,
    }
    return rows.reshape(-1, 6), summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--modes", nargs="+", default=None, help="default: every mode in the manifest")
    ap.add_argument("--max-events", type=int, default=None)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--work-dir", type=Path, default=None)
    args = ap.parse_args(argv)

    work = args.work_dir or mg5_work_dir()
    manifest = json.loads((work / "higgs_lhe_manifest.json").read_text())
    modes = args.modes or list(manifest)
    out_manifest_path = work / "higgs_shower_manifest.json"
    out_manifest = json.loads(out_manifest_path.read_text()) if out_manifest_path.exists() else {}
    for mode in modes:
        lhe = Path(manifest[mode]["lhe"])
        print(f"[{mode}] showering {lhe.name} ...", flush=True)
        rows, summary = shower_lhe(lhe, seed=args.seed, max_events=args.max_events)
        csv = work / f"higgs_{mode}.csv"
        np.savetxt(csv, rows, delimiter=",", fmt="%.8e", header="pt_GeV,y,E,px,py,pz", comments="# ")
        summary.update({"csv": str(csv), "lhe": str(lhe), "sigma_mg5_lo_pb": manifest[mode].get("sigma_lo_pb"),
                        "timestamp": datetime.now().isoformat()})
        out_manifest[mode] = summary
        out_manifest_path.write_text(json.dumps(out_manifest, indent=2) + "\n")
        print(f"[{mode}] {summary['n_higgs']} Higgs, <pT> = {rows[:, 0].mean():.1f} GeV, "
              f"<|y|> = {np.abs(rows[:, 1]).mean():.2f}, sigma(LHE) = {summary['sigma_lhe_pb']:.3f} pb", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
