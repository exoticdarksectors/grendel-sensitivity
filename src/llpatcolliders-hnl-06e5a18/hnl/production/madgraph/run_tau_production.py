#!/usr/bin/env python3
"""Prompt-tau production: MG5 tau pool, then tau decays to HNL CSVs."""

import argparse
import random
import shutil
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config_mass_grid import MASS_GRID, N_EVENTS_DEFAULT
from production.constants import K_FACTOR_EW_BY_PROCESS, M_TAU
from production.decay_engine.tau_decay import (
    init_hnlcalc, compute_tau_production_br_components, sample_hnl_from_tau,
    TAU_W_2BODY_ASYMMETRY,
)
from production.madgraph._mg5_common import MG5_EXE, LHAPDF_CONFIG
from production.madgraph.runner import (
    ensure_process_dir, run_events as mg5_run_events, write_run_card,
)
from production.madgraph.lhe_to_csv import LHEParser
from production.io import (
    llp_csv_path, read_csv_matrix, write_csv_matrix, write_empty_csv,
)
from production.paths import (
    MG5_WORK_DIR, TAU_POOL_CSV, existing_tau_pool_csv,
)

CARDS_DIR = Path(__file__).parent / "cards"
WORK_DIR = MG5_WORK_DIR

POOL_CSV = TAU_POOL_CSV


def _build_process_dir():
    return ensure_process_dir(
        label="tau_pool",
        model_import="import model sm",
        proc_card=CARDS_DIR / "proc_card_tau_production.dat",
        work_dir=WORK_DIR,
        process_dir=WORK_DIR / "tau_pool",
        generation_timeout=600,
    )


def generate_tau_pool(n_events, nb_core=1, pool_path=None):
    if pool_path is None:
        pool_path = POOL_CSV
    print(f"\nStage 1: generating tau pool ({n_events} events, {nb_core} core(s)) → {pool_path}")

    work_subdir = _build_process_dir()
    if work_subdir is None:
        return False

    write_run_card(work_subdir, CARDS_DIR, n_events)

    lhe_path = mg5_run_events(work_subdir, "run_pool", nb_core=nb_core, timeout=7200)
    if lhe_path is None:
        return False

    pool_path.parent.mkdir(parents=True, exist_ok=True)
    n_tau = LHEParser(lhe_path).write_tau_csv(pool_path)
    if n_tau == 0:
        print("    FAILED: no taus extracted from LHE")
        return False
    print(f"    OK: {n_tau} tau rows → {pool_path}")

    # The pool carries the parent origin in column 6 (w, E, px, py, pz, origin):
    # |origin| == 24 marks W-origin taus, everything else is Drell-Yan
    # (gamma*/Z -> tau tau). Scale each row by its process K-factor.
    data = read_csv_matrix(pool_path)
    if data.size:
        is_w = np.abs(data[:, 5]) == 24
        data[is_w, 0] *= K_FACTOR_EW_BY_PROCESS["W"]
        data[~is_w, 0] *= K_FACTOR_EW_BY_PROCESS["DY"]
        write_csv_matrix(pool_path, data)

    shutil.rmtree(work_subdir / "Events" / "run_pool", ignore_errors=True)
    return True


def process_flavor(flavor, pool_w, pool_E, pool_px, pool_py, pool_pz, pool_origin,
                   masses, rng):
    hnl = init_hnlcalc(flavor)

    # W-origin taus are polarized (natural helicity: the 2-body N comes out
    # forward, asym = +1); all DY-like rows are treated unpolarized.
    is_W = np.abs(pool_origin) == 24
    is_DY = ~is_W

    for m_N in masses:
        csv_path = llp_csv_path(flavor, "tau", m_N)

        if m_N >= M_TAU:
            write_empty_csv(csv_path)
            continue

        br_2body, br_3body, br_total = compute_tau_production_br_components(hnl, m_N)
        if br_total <= 0:
            write_empty_csv(csv_path)
            continue

        out_rows = []
        for sel, asym in (
            (is_W,  TAU_W_2BODY_ASYMMETRY),
            (is_DY, 0.0),
        ):
            if not sel.any():
                continue
            hnl_4v, _, _ = sample_hnl_from_tau(
                pool_E[sel], pool_px[sel], pool_py[sel], pool_pz[sel],
                m_N, br_2body, br_3body, br_total, rng, asymmetry=asym,
            )
            weights = pool_w[sel] * br_total
            block = np.column_stack([
                weights, hnl_4v[:, 0], hnl_4v[:, 1], hnl_4v[:, 2], hnl_4v[:, 3],
            ])
            out_rows.append(block)

        if out_rows:
            data = np.vstack(out_rows)
            write_csv_matrix(csv_path, data)
            print(f"    {csv_path.name}: {len(data)} events, BR_total={br_total:.3e}, "
                  f"w_sum={data[:, 0].sum():.3e}")
        else:
            write_empty_csv(csv_path)


def _load_pool(pool_path=None):
    if pool_path is None:
        pool_path = existing_tau_pool_csv()
    if not pool_path.exists() or pool_path.stat().st_size == 0:
        raise FileNotFoundError(
            f"Tau pool not found at {pool_path}. Run without --skip-mg5 (Stage 1) "
            f"to generate it."
        )
    data = read_csv_matrix(pool_path)
    if data.shape[1] != 6:
        raise ValueError(
            f"Tau pool {pool_path} has {data.shape[1]} columns; expected 6 "
            f"(w, E, px, py, pz, origin). Regenerate with the current driver."
        )
    return data.T


def main():
    ap = argparse.ArgumentParser(description="Prompt-tau HNL production (W → τν, Z → ττ)")
    ap.add_argument("--flavor", choices=["Ue", "Umu", "Utau"], nargs="+",
                    default=["Ue", "Umu", "Utau"])
    ap.add_argument("--masses", type=float, nargs="+", default=None,
                    help="Custom mass list (default: full grid, m < m_tau)")
    ap.add_argument("--nevents", type=int, default=N_EVENTS_DEFAULT,
                    help="Stage 1 tau pool size (default: %(default)d)")
    ap.add_argument("--nb-core", type=int, default=1,
                    help="MG5 cores for Stage 1")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--skip-mg5", action="store_true",
                    help="Skip Stage 1; reuse cached tau_pool.csv")
    ap.add_argument("--test", action="store_true",
                    help="Smoke test: Umu only, single mass (1 GeV), 1k pool")
    args = ap.parse_args()

    if args.test:
        flavors = ["Umu"]
        masses = [1.0]
        n_events = 1000
        pool_path = WORK_DIR / "tau_pool_test.csv"
    else:
        flavors = args.flavor
        masses = args.masses if args.masses else MASS_GRID
        n_events = args.nevents
        pool_path = existing_tau_pool_csv() if args.skip_mg5 else POOL_CSV

    rng = np.random.default_rng(args.seed)
    random.seed(args.seed)

    if not MG5_EXE.exists() and not args.skip_mg5:
        print(f"ERROR: MadGraph not found at {MG5_EXE}")
        print("Set $HNL_MG5_EXE or vendor MG5 under hnl/vendored/ "
              "(see hnl/vendored/PROVENANCE.md).")
        return 1
    if not LHAPDF_CONFIG.exists() and not args.skip_mg5:
        print(f"ERROR: lhapdf-config not found at {LHAPDF_CONFIG}")
        return 1

    print(f"Prompt-tau production")
    print(f"  MG5: {MG5_EXE}")
    print(f"  LHAPDF: {LHAPDF_CONFIG}")
    print(f"  Pool CSV: {pool_path}")
    print(f"  Flavors: {flavors}")
    print(f"  Masses: {len(masses)} points")
    print(f"  Stage 1: {'skipped (using cached pool)' if args.skip_mg5 else f'{n_events} events'}")

    if not args.skip_mg5:
        if not generate_tau_pool(n_events, nb_core=args.nb_core, pool_path=pool_path):
            print("\nStage 1 failed. Aborting.")
            return 1

    pool_w, pool_E, pool_px, pool_py, pool_pz, pool_origin = _load_pool(pool_path)
    n_W  = int((np.abs(pool_origin) == 24).sum())
    n_DY = len(pool_origin) - n_W
    print(f"\nStage 2: decaying pool ({len(pool_origin)} taus: W±={n_W} polarised, "
          f"γ*/Z DY={n_DY} unpolarised)")

    for flavor in flavors:
        print(f"\nFlavor: {flavor}")
        process_flavor(flavor, pool_w, pool_E, pool_px, pool_py, pool_pz, pool_origin,
                       masses, rng)

    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
