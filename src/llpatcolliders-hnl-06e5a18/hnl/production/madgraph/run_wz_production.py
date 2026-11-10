#!/usr/bin/env python3
"""MadGraph W/Z -> l N production for HNL CSVs."""

import sys
import shutil
import argparse
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config_mass_grid import MASS_GRID, N_EVENTS_DEFAULT
from production.constants import K_FACTOR_EW_BY_PROCESS, FLAVOR_TO_MG5
from production.madgraph._mg5_common import MG5_EXE, LHAPDF_CONFIG
from production.madgraph.runner import (
    ensure_process_dir, run_events as mg5_run_events, write_run_card,
)
from production.io import llp_csv_path, scale_weight_column
from production.paths import MG5_WORK_DIR

MODEL_DIR = PROJECT_ROOT / "vendored" / "SM_HeavyN_CKM_AllMasses_LO"

CARDS_DIR = Path(__file__).parent / "cards"
WORK_DIR = MG5_WORK_DIR

MIXING_CONFIGS = {
    "Ue":   {"ve1": 1.0, "vmu1": 0.0, "vtau1": 0.0},
    "Umu":  {"ve1": 0.0, "vmu1": 1.0, "vtau1": 0.0},
    "Utau": {"ve1": 0.0, "vmu1": 0.0, "vtau1": 1.0},
}


def get_or_create_process_dir(flavor):
    mg5_flavor = FLAVOR_TO_MG5[flavor]
    proc_card = CARDS_DIR / f"proc_card_{mg5_flavor}.dat"
    return ensure_process_dir(
        label=mg5_flavor,
        model_import=f"import model {MODEL_DIR}",
        proc_card=proc_card,
        work_dir=WORK_DIR,
        process_dir=WORK_DIR / f"hnl_{mg5_flavor}",
        generation_timeout=300,
    )


def write_cards(work_subdir, flavor, mass, n_events):
    cards_dir = work_subdir / "Cards"
    cards_dir.mkdir(exist_ok=True)

    write_run_card(work_subdir, CARDS_DIR, n_events)

    param_content = (CARDS_DIR / "param_card_template.dat").read_text()
    mixing = MIXING_CONFIGS[flavor]
    param_content = param_content.replace("MASS_N1_PLACEHOLDER", f"{mass:.6e}")
    param_content = param_content.replace("VE1_PLACEHOLDER", f"{mixing['ve1']:.6e}")
    param_content = param_content.replace("VMU1_PLACEHOLDER", f"{mixing['vmu1']:.6e}")
    param_content = param_content.replace("VTAU1_PLACEHOLDER", f"{mixing['vtau1']:.6e}")
    (cards_dir / "param_card.dat").write_text(param_content)


def convert_lhe(lhe_path, csv_path):
    from production.madgraph.lhe_to_csv import LHEParser
    parser = LHEParser(lhe_path)
    return parser.write_hnl_csv(csv_path)


def run_single_point(flavor, mass, n_events, nb_core=1):
    csv_path = llp_csv_path(flavor, "WZ", mass)
    csv_path.unlink(missing_ok=True)

    print(f"\n  [{flavor}] m_N = {mass} GeV")

    work_subdir = get_or_create_process_dir(flavor)
    if work_subdir is None:
        return False

    write_cards(work_subdir, flavor, mass, n_events)

    run_name = f"run_{csv_path.stem.removeprefix('mN_')}"
    lhe_path = mg5_run_events(work_subdir, run_name, nb_core=nb_core)
    if lhe_path is None:
        return False

    n_ev = convert_lhe(lhe_path, csv_path)
    if n_ev is None or n_ev == 0:
        print(f"    FAILED: no HNL events extracted")
        return False

    print(f"    OK: {n_ev} events → {csv_path}")

    # The sample mixes p p > W -> l N (dominant) with p p > Z -> nu N in one
    # LHE/CSV and carries no per-row process tag, so the whole sample is scaled
    # by the W K-factor. A per-row W vs Z split would require carrying the LHE
    # mother PDG into the HNL rows (a future refinement; see
    # `REMAINING_WORK.md`).
    scale_weight_column(csv_path, K_FACTOR_EW_BY_PROCESS["W"])

    shutil.rmtree(work_subdir / "Events" / run_name, ignore_errors=True)
    return True


def main():
    parser = argparse.ArgumentParser(description="W/Z → ℓ N production via MadGraph")
    parser.add_argument("--flavor", choices=["Ue", "Umu", "Utau"], nargs="+",
                        default=["Ue", "Umu", "Utau"])
    parser.add_argument("--masses", type=float, nargs="+", default=None)
    parser.add_argument("--nevents", type=int, default=N_EVENTS_DEFAULT)
    parser.add_argument("--nb-core", type=int, default=1,
                        help="Number of CPU cores per MG5 generate_events job")
    parser.add_argument("--test", action="store_true",
                        help="Quick test: single point (Umu, 1.0 GeV, 1000 events)")
    parser.add_argument("--min-mass", type=float, default=None)
    args = parser.parse_args()

    if args.test:
        flavors = ["Umu"]
        masses = [1.0]
        n_events = 1000
    else:
        flavors = args.flavor
        masses = args.masses if args.masses else MASS_GRID
        n_events = args.nevents

    nb_core = args.nb_core
    if nb_core < 1:
        print("ERROR: --nb-core must be >= 1")
        return 1

    if args.min_mass is not None:
        masses = [m for m in masses if m >= args.min_mass]

    if not MG5_EXE.exists():
        print(f"ERROR: MadGraph not found at {MG5_EXE}")
        print("Set $HNL_MG5_EXE or vendor MG5 under hnl/vendored/ "
              "(see hnl/vendored/PROVENANCE.md).")
        return 1
    if not MODEL_DIR.exists():
        print(f"ERROR: HeavyN UFO model not found at {MODEL_DIR}")
        return 1
    if not LHAPDF_CONFIG.exists():
        print(f"ERROR: lhapdf-config not found at {LHAPDF_CONFIG}")
        print("Set $HNL_LHAPDF_CONFIG to a working lhapdf-config binary "
              "(the LHAPDF install must carry NNPDF40_nlo_as_01180).")
        return 1

    print(f"W/Z → ℓ N Production")
    print(f"  MG5: {MG5_EXE}")
    print(f"  Model: {MODEL_DIR}")
    print(f"  LHAPDF: {LHAPDF_CONFIG}")
    print(f"  Flavors: {flavors}")
    print(f"  Masses: {len(masses)} points")
    print(f"  Events/point: {n_events}")
    print(f"  Cores/point: {nb_core}")

    n_ok = 0
    n_fail = 0
    for flavor in flavors:
        for mass in masses:
            if run_single_point(flavor, mass, n_events, nb_core=nb_core):
                n_ok += 1
            else:
                n_fail += 1

    print(f"\nDone: {n_ok} OK, {n_fail} failed")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
