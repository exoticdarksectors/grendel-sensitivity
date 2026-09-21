#!/usr/bin/env python3
"""Do the exHad-decayed curves fall inside the variation bands?

The uncertainty bundles are produced with the primary decay models: the
FairShip exclusive-channel table (HNL), the analytic two-body Winkler proxy
(BC4) and the Pythia templates (BC10). Each spans its decay-model axis
differently; the exHad model differs from every one of them in both width
table and final states, so the first question is whether a band already
covers the exHad curve. Where it does not, the supplements
(``grendel.band.hnl.decay_model_band_exhad``, ``grendel.band.scalar.exhad_variation``)
add exHad to the axis the way each band was built.

Both edges are compared with the band ordered (``min(lo, hi)`` to
``max(lo, hi)``): the HNL ``*_dm_lo/hi`` columns are named after the width
direction, not the vertical order, and cross above 3 GeV.

    python scripts/check_curves_against_bands.py --results-dir RESULTS

expects ``RESULTS/<model>/sensitivity_exhad.csv`` next to each model's
``sensitivity.csv`` and the band CSVs under ``RESULTS/<model>/band/``.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def col(frame: pd.DataFrame, name: str):
    """The exHad edge column, whether or not the merge had to suffix it."""
    return frame[f"{name}_ex"] if f"{name}_ex" in frame else frame[name]


def report(title: str, value, lo, hi) -> dict:
    value = np.asarray(value, float)
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    ok = np.isfinite(value) & np.isfinite(lo) & np.isfinite(hi)
    value = value[ok]
    lo, hi = np.fmin(lo[ok], hi[ok]), np.fmax(lo[ok], hi[ok])
    inside = (value >= lo) & (value <= hi)
    # Signed distance past whichever edge was crossed, in decades.
    excess = np.where(value > hi, np.log10(value / hi), np.log10(value / lo))
    excess = np.where(inside, 0.0, excess)
    worst = float(np.abs(excess).max()) if len(excess) else float("nan")
    print(f"  {title:<44} {int(inside.sum()):>3}/{len(inside):<3} inside "
          f"({100 * inside.mean() if len(inside) else 0:3.0f}%)   worst breach {worst:.3f} dex")
    return {"n": int(len(inside)), "inside": int(inside.sum()), "worst_dex": worst}


def sensitive(merged: pd.DataFrame, a: str, b: str) -> pd.DataFrame:
    return merged[merged[a].astype(bool) & merged[b].astype(bool)]


def check_hnl(results: Path) -> None:
    curve = results / "hnl" / "sensitivity_exhad.csv"
    band = results / "hnl" / "band" / "decay_model_band.csv"
    if not (curve.exists() and band.exists()):
        print("\nHNL: skipped (no exHad curve or decay_model_band.csv)")
        return
    x = pd.read_csv(curve)
    print("\nHNL -- decay_model_band.csv (+/-delta width legs)")
    b = pd.read_csv(band)
    m = sensitive(b.merge(x, on=["mass_GeV", "flavor"], suffixes=("_b", "_ex")), "has_sensitivity_b", "has_sensitivity_ex")
    report("lower edge", col(m, "u2_min"), m.u2_min_dm_lo, m.u2_min_dm_hi)
    report("upper edge", col(m, "u2_max"), m.u2_max_dm_lo, m.u2_max_dm_hi)
    supplement = results / "hnl" / "band" / "decay_model_band_exhad.csv"
    if supplement.exists():
        print("HNL -- decay_model_band_exhad.csv (legs + exHad member)")
        s = pd.read_csv(supplement)
        m = sensitive(s.merge(x, on=["mass_GeV", "flavor"], suffixes=("_b", "_ex")), "has_sensitivity_b", "has_sensitivity_ex")
        report("lower edge", col(m, "u2_min"), m.u2_min_dm_lo, m.u2_min_dm_hi)
        report("upper edge", col(m, "u2_max"), m.u2_max_dm_lo, m.u2_max_dm_hi)


def check_bc4(results: Path) -> None:
    curve = results / "bc4" / "sensitivity_exhad.csv"
    band = results / "bc4" / "band" / "bc4_single_source_variation_envelope.csv"
    if not (curve.exists() and band.exists()):
        print("\nBC4: skipped (no exHad curve or envelope)")
        return
    x = pd.read_csv(curve)
    for label, path in (("campaign", band), ("supplemented", band.with_name("bc4_single_source_variation_envelope_exhad.csv"))):
        if not path.exists():
            continue
        print(f"\nBC4 -- {path.name} ({label})")
        e = pd.read_csv(path)
        m = sensitive(e.merge(x, on="mass_GeV", suffixes=("_env", "_ex")), "has_sensitivity_env", "has_sensitivity_ex")
        report("lower edge, full envelope", col(m, "u2_min"), m.u2_min_envelope_lo, m.u2_min_envelope_hi)
        report("lower edge, decay axis only", col(m, "u2_min"), m.u2_min_decay_model_envelope_lo, m.u2_min_decay_model_envelope_hi)
        report("upper edge, full envelope", col(m, "u2_max"), m.u2_max_envelope_lo, m.u2_max_envelope_hi)
        report("upper edge, decay axis only", col(m, "u2_max"), m.u2_max_decay_model_envelope_lo, m.u2_max_decay_model_envelope_hi)


def check_bc10(results: Path) -> None:
    curve = results / "bc10" / "sensitivity_exhad.csv"
    band = results / "bc10" / "band" / "bc10_single_source_variation_envelope.csv"
    if not (curve.exists() and band.exists()):
        print("\nBC10: skipped (no exHad curve or envelope)")
        return
    print("\nBC10 -- bc10_single_source_variation_envelope.csv")
    e, x = pd.read_csv(band), pd.read_csv(curve)
    m = sensitive(e.merge(x, on="mass_GeV", suffixes=("_env", "_ex")), "has_sensitivity_env", "has_sensitivity_ex")
    report("lower edge, full envelope", col(m, "invf_min"), m.invf_min_envelope_lo, m.invf_min_envelope_hi)
    report("upper edge, full envelope", col(m, "invf_max"), m.invf_max_envelope_lo, m.invf_max_envelope_hi)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", type=Path, required=True)
    args = ap.parse_args()
    print(__doc__.splitlines()[0])
    check_hnl(args.results_dir)
    check_bc4(args.results_dir)
    check_bc10(args.results_dir)
    print("""
Reading these numbers: the exHad curves are re-solved on the same parent pools
with only the decay stage changed, so they are the exHad model at the same
production. A supplemented band contains that curve up to the Monte-Carlo
scatter between the supplement's own run and the campaign's, so residual
breaches at the 0.01-0.05 dex level are numerical, not model, effects.""")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
