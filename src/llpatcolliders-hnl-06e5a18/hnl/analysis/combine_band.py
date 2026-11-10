#!/usr/bin/env python3
"""Combine per-variation HNL sensitivity curves into a theory-uncertainty band.

Reads the ``band_registry.json`` written by ``run_variation_band.py`` (one
``hnl_sensitivity.csv`` per coherent FONLL variation) and turns the spread of
each exclusion boundary into an uncertainty ribbon, per the agreed prescription:

* all combination is done in ``x = log10(U^2)`` (the boundary spans many
  decades, so min/max or std in linear U^2 would be meaningless);
* **scale** -> asymmetric envelope of the 7 coherent scale curves
  (max/min deviation of ``x`` from central), per boundary;
* **PDF** -> ``std(x, ddof=1)`` over the replica members (the NNPDF
  Monte-Carlo 1 sigma), per boundary;
* **mass** -> bottom and charm are independent: the per-quark maximum
  ``|x - x_central|`` are added in quadrature;
* the three axes are combined in quadrature (uncorrelated) into an up/down
  ``sigma_x`` per boundary; the ribbon is ``10**(x_central +/- sigma_x)``.

Open-edge (``u2_*_open``) and no-sensitivity states are preserved: a boundary
is reported open if central or any contributing variation ran off the scan
edge, and masses where central has no sensitivity carry no band.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HNL_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_REGISTRY = HNL_ROOT / "tmp" / "runs" / "band_registry.json"

BOUNDARIES = [("u2_min", "u2_min_open"), ("u2_max", "u2_max_open")]


def _load_curves(registry: dict) -> dict[str, pd.DataFrame]:
    curves = {}
    for v in registry["variations"]:
        df = pd.read_csv(v["sensitivity_csv"])
        df = df.set_index(["flavor", "mass_GeV"]).sort_index()
        curves[v["name"]] = df
    return curves


def _axis_members(registry: dict) -> dict[str, list[str]]:
    axes: dict[str, list[str]] = {"central": [], "scale": [], "pdf": [], "mass": []}
    for v in registry["variations"]:
        axes[v["axis"]].append(v["name"])
    return axes


def _boundary(curve: pd.DataFrame, key, mass_key, col, open_col):
    """Return (x=log10(boundary), is_open, has_sens) for one variation/point."""
    if key not in curve.index:
        return None, False, False
    row = curve.loc[key]
    if not bool(row.get("has_sensitivity", False)):
        return None, False, False
    is_open = bool(row.get(open_col, False))
    val = row.get(col, np.nan)
    if not np.isfinite(val) or val <= 0:
        return None, is_open, True
    return float(np.log10(val)), is_open, True


def combine_band(registry_path: Path) -> pd.DataFrame:
    registry = json.loads(registry_path.read_text())
    curves = _load_curves(registry)
    axes = _axis_members(registry)
    if not axes["central"]:
        raise ValueError("registry has no central variation")
    central_name = axes["central"][0]
    central = curves[central_name]

    out_rows = []
    for (flavor, mass) in central.index:
        crow = central.loc[(flavor, mass)]
        rec = {"flavor": flavor, "mass_GeV": mass,
               "has_sensitivity": bool(crow.get("has_sensitivity", False))}
        if not rec["has_sensitivity"]:
            out_rows.append(rec)
            continue

        for col, open_col in BOUNDARIES:
            xc, c_open, _ = _boundary(central, (flavor, mass), mass, col, open_col)
            rec[f"{col}_central"] = crow.get(col, np.nan)
            rec[f"{col}_open"] = c_open
            if xc is None:
                # central boundary open/undefined -> ribbon open on this side
                rec[f"{col}_band_lo"] = np.nan
                rec[f"{col}_band_hi"] = np.nan
                rec[f"{col}_open"] = True
                continue

            any_open = c_open

            # --- scale: asymmetric envelope of central + scale curves ---
            xs = [xc]
            for name in axes["scale"]:
                x, op, _ = _boundary(curves[name], (flavor, mass), mass, col, open_col)
                if op:
                    any_open = True
                if x is not None:
                    xs.append(x)
            scale_up = max(xs) - xc
            scale_dn = xc - min(xs)

            # --- pdf: std over replica members in log10(U^2) ---
            xp = []
            for name in axes["pdf"]:
                x, op, _ = _boundary(curves[name], (flavor, mass), mass, col, open_col)
                if op:
                    any_open = True
                if x is not None:
                    xp.append(x)
            pdf_sigma = float(np.std(xp, ddof=1)) if len(xp) >= 2 else 0.0

            # --- mass: independent bottom/charm, per-quark max |dev|, in quadrature ---
            mdev = {"mb": 0.0, "mc": 0.0}
            for name in axes["mass"]:
                x, op, _ = _boundary(curves[name], (flavor, mass), mass, col, open_col)
                if op:
                    any_open = True
                if x is None:
                    continue
                quark = "mb" if name.startswith("mb") else "mc"
                mdev[quark] = max(mdev[quark], abs(x - xc))
            mass_dev = float(np.hypot(mdev["mb"], mdev["mc"]))

            sigma_up = float(np.sqrt(scale_up**2 + pdf_sigma**2 + mass_dev**2))
            sigma_dn = float(np.sqrt(scale_dn**2 + pdf_sigma**2 + mass_dev**2))

            rec[f"{col}_band_lo"] = 10.0 ** (xc - sigma_dn)
            rec[f"{col}_band_hi"] = 10.0 ** (xc + sigma_up)
            rec[f"{col}_open"] = any_open
            rec[f"{col}_scale_up_dex"] = scale_up
            rec[f"{col}_scale_dn_dex"] = scale_dn
            rec[f"{col}_pdf_sigma_dex"] = pdf_sigma
            rec[f"{col}_mass_dev_dex"] = mass_dev

        out_rows.append(rec)

    return pd.DataFrame(out_rows).sort_values(["flavor", "mass_GeV"])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--registry", default=str(DEFAULT_REGISTRY),
                    help="band_registry.json from run_variation_band.py")
    ap.add_argument("--out", default=None,
                    help="output band CSV (default: <registry dir>/hnl_band.csv)")
    args = ap.parse_args(argv)

    registry_path = Path(args.registry)
    if not registry_path.exists():
        print(f"registry not found: {registry_path}; run run_variation_band.py first")
        return 1

    df = combine_band(registry_path)
    out = Path(args.out) if args.out else registry_path.parent / "hnl_band.csv"
    df.to_csv(out, index=False)

    reg = json.loads(registry_path.read_text())
    axes = _axis_members(reg)
    n_sens = int(df["has_sensitivity"].sum())
    print(f"combined {len(reg['variations'])} variations "
          f"(scale={len(axes['scale'])}, pdf={len(axes['pdf'])}, mass={len(axes['mass'])})")
    print(f"band written: {out}  ({n_sens}/{len(df)} mass points with sensitivity)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
