"""BC5 Higgs production, stage 1: the hard process of every SM mode with MadGraph5_aMC@NLO, one LHE
file per mode."""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from datetime import datetime
from pathlib import Path

from ...io.paths import work_dir
from ..hnl.production.madgraph._mg5_common import LHAPDF_CONFIG, MG5_EXE
from ..hnl.production.madgraph.runner import ensure_process_dir, run_events, write_run_card

CARDS_DIR = Path(__file__).resolve().parent / "cards"
MODES = {
    "ggf": "import model heft",
    "vbf": "import model sm",
    "wh": "import model sm",
    "zh": "import model sm",
    "tth": "import model sm",
}
N_EVENTS_DEFAULT = 200_000


def mg5_work_dir() -> Path:
    value = os.environ.get("GRENDEL_BC5_MG5_WORK_DIR")
    return Path(value).expanduser() if value else work_dir() / "bc5" / "madgraph"


def lhe_cross_section_pb(lhe_path: Path) -> float | None:
    """The MG5 LO cross section written in the LHE banner."""
    import gzip
    opener = gzip.open if lhe_path.suffix == ".gz" else open
    with opener(lhe_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            match = re.search(r"Integrated weight \(pb\)\s*:\s*([0-9.eE+-]+)", line)
            if match:
                return float(match.group(1))
            if line.startswith("<event"):
                break
    return None


def run_mode(mode: str, n_events: int, nb_core: int, work: Path) -> dict | None:
    process_dir = ensure_process_dir(
        label=f"higgs_{mode}", model_import=MODES[mode],
        proc_card=CARDS_DIR / f"proc_card_{mode}.dat", work_dir=work,
        process_dir=work / f"higgs_{mode}", generation_timeout=1800, compile_timeout=1800)
    if process_dir is None:
        return None
    write_run_card(process_dir, CARDS_DIR, n_events)
    run_name = f"run_{n_events}"
    lhe = run_events(process_dir, run_name, nb_core=nb_core, timeout=4 * 3600)
    if lhe is None:
        return None
    return {"mode": mode, "lhe": str(lhe), "n_events_requested": n_events,
            "sigma_lo_pb": lhe_cross_section_pb(lhe)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--modes", nargs="+", choices=sorted(MODES), default=list(MODES))
    ap.add_argument("--nevents", type=int, default=N_EVENTS_DEFAULT)
    ap.add_argument("--nb-core", type=int, default=1)
    ap.add_argument("--work-dir", type=Path, default=None,
                    help="MG5 work directory (default: $GRENDEL_BC5_MG5_WORK_DIR)")
    args = ap.parse_args(argv)

    if not MG5_EXE.exists():
        print(f"ERROR: MadGraph not found at {MG5_EXE}; set $GRENDEL_MG5_EXE")
        return 1
    if not LHAPDF_CONFIG.exists():
        print(f"ERROR: lhapdf-config not found at {LHAPDF_CONFIG}; set $GRENDEL_LHAPDF_CONFIG")
        return 1
    work = args.work_dir or mg5_work_dir()
    work.mkdir(parents=True, exist_ok=True)
    manifest_path = work / "higgs_lhe_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    print(f"MG5: {MG5_EXE}\nLHAPDF: {LHAPDF_CONFIG}\nwork: {work}")
    failed = []
    for mode in args.modes:
        print(f"\n[{mode}] {args.nevents} events ...", flush=True)
        entry = run_mode(mode, args.nevents, args.nb_core, work)
        if entry is None:
            failed.append(mode)
            continue
        entry["timestamp"] = datetime.now().isoformat()
        manifest[mode] = entry
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        print(f"[{mode}] sigma_LO = {entry['sigma_lo_pb']} pb -> {entry['lhe']}", flush=True)
    if failed:
        print(f"\nFAILED modes: {failed}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
