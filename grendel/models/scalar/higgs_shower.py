"""BC5 Higgs production, stage 2: shower the MadGraph LHE events with Pythia 8 and record the Higgs
four-vector after the shower."""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys
import tempfile
from datetime import datetime
from pathlib import Path

import numpy as np

from .higgs_production import mg5_work_dir

HIGGS_PDG = 25


def _plain_lhe(lhe: Path, tmp_dir: Path) -> Path:
    """Pythia's LHEF reader wants a plain file unless built with gzip."""
    if lhe.suffix != ".gz":
        return lhe
    out = tmp_dir / lhe.with_suffix("").name
    with gzip.open(lhe, "rb") as src, open(out, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return out


def shower_lhe(lhe: Path, *, seed: int = 1, max_events: int | None = None, quiet: bool = True):
    import pythia8

    with tempfile.TemporaryDirectory() as tmp:
        plain = _plain_lhe(lhe, Path(tmp))
        pythia = pythia8.Pythia("", False)
        for line in (
            "Beams:frameType = 4",
            f"Beams:LHEF = {plain}",
            "PartonLevel:ISR = on",
            "PartonLevel:FSR = on",
            "PartonLevel:MPI = off",
            "HadronLevel:all = off",
            f"{HIGGS_PDG}:mayDecay = off",
            "ProcessLevel:resonanceDecays = off",
            "Check:event = off",
            "Next:numberCount = 0",
            "Next:numberShowEvent = 0",
            "Next:numberShowInfo = 0",
            "Next:numberShowProcess = 0",
            "Next:numberShowLHA = 0",
            "Init:showChangedSettings = off",
            "Init:showChangedParticleData = off",
            f"Print:quiet = {'on' if quiet else 'off'}",
            "Random:setSeed = on",
            f"Random:seed = {seed}",
        ):
            pythia.readString(line)
        if not pythia.init():
            raise RuntimeError(f"Pythia init failed on {lhe}")
        rows = []
        n_done = 0
        while max_events is None or n_done < max_events:
            if not pythia.next():
                if pythia.infoPython().atEndOfFile():
                    break
                continue
            n_done += 1
            higgs = None
            for particle in pythia.event:
                if particle.id() == HIGGS_PDG and particle.isFinal():
                    higgs = particle
            if higgs is None:
                continue
            rows.append((higgs.pT(), higgs.y(), higgs.e(), higgs.px(), higgs.py(), higgs.pz()))
        info = pythia.infoPython()
        summary = {
            "pythia_version": pythia.parm("Pythia:versionNumber"),
            "n_showered": n_done,
            "n_higgs": len(rows),
            "sigma_lhe_pb": info.sigmaGen() * 1e9,
            "seed": seed,
        }
    return np.array(rows, float).reshape(-1, 6), summary


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
