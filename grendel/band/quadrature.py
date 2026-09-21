"""The HNL production band: FONLL variations combined in quadrature.

One sensitivity curve per coherent FONLL variation (the same scale choice or
PDF replica applied to bottom and charm together; the two quark masses
varied independently) is turned into an uncertainty ribbon on each
exclusion boundary, in ``x = log10(U^2)``:

* **scale** -- asymmetric envelope of the central and the six scale curves
  (largest deviation of ``x`` up and down);
* **PDF** -- ``std(x, ddof=1)`` over the replica members, the NNPDF
  Monte-Carlo one-sigma;
* **mass** -- bottom and charm are independent: the per-quark maximum
  ``|x - x_central|`` added in quadrature;
* **alpha_s** (optional, PDF4LHC) -- half the ``|x|`` spread between the
  as = 0.119 and 0.117 companion curves, folded symmetrically;
* the axes are combined in quadrature into ``sigma_up``/``sigma_dn`` per
  boundary; the ribbon is ``10**(x_central +/- sigma)``.

Open edges and no-sensitivity states are preserved: a boundary is open if
the central or any member ran off the scan edge, and masses where the
central has no sensitivity carry no band. If any member creates or
destroys an island relative to the central, the point is marked
``*_topology_changed`` and no ribbon is reported -- a topology change is
not reduced to a band from the surviving members.

``suppress`` names ``(flavor, mass)`` points whose ``u2_min`` ribbon is
withheld after combination (``band_lo/hi`` NaN, ``band_suppressed`` True)
because the member spread there is acceptance-Monte-Carlo noise rather
than a resolvable production uncertainty; the central and the dex
components stay for auditability.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

BOUNDARIES = [("u2_min", "u2_min_open"), ("u2_max", "u2_max_open")]
AXES = ("central", "scale", "pdf", "mass")


def _boundary(curve: pd.DataFrame, key, col, open_col):
    """(x = log10 boundary, is_open, has_sensitivity) for one member at one point."""
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


def read_curve(path) -> pd.DataFrame:
    """A sensitivity CSV indexed by (flavor, mass_GeV)."""
    return pd.read_csv(path).set_index(["flavor", "mass_GeV"]).sort_index()


def combine_curves(curves: dict[str, pd.DataFrame], axes: dict[str, list[str]], *,
                   alphas: tuple[pd.DataFrame, pd.DataFrame] | None = None,
                   suppress=frozenset()) -> pd.DataFrame:
    """Combine indexed member curves (see :func:`read_curve`) into the ribbon.

    ``axes`` maps ``central``/``scale``/``pdf``/``mass`` to member names;
    ``alphas`` is the ``(as_lo, as_hi)`` companion pair or None.
    """
    if not axes.get("central"):
        raise ValueError("no central variation")
    central = curves[axes["central"][0]]
    members = list(axes.get("scale", [])) + list(axes.get("pdf", [])) + list(axes.get("mass", []))
    as_lo, as_hi = alphas if alphas is not None else (None, None)
    fold_alphas = as_lo is not None and as_hi is not None
    suppress = {(flavor, float(mass)) for flavor, mass in suppress}

    out_rows = []
    for (flavor, mass) in central.index:
        crow = central.loc[(flavor, mass)]
        rec = {"flavor": flavor, "mass_GeV": mass,
               "has_sensitivity": bool(crow.get("has_sensitivity", False))}
        if not rec["has_sensitivity"]:
            out_rows.append(rec)
            continue

        for col, open_col in BOUNDARIES:
            xc, c_open, _ = _boundary(central, (flavor, mass), col, open_col)
            rec[f"{col}_central"] = crow.get(col, np.nan)
            rec[f"{col}_open"] = c_open
            member_states = [_boundary(curves[name], (flavor, mass), col, open_col)
                             for name in members]
            if fold_alphas:
                member_states.extend([_boundary(as_lo, (flavor, mass), col, open_col),
                                      _boundary(as_hi, (flavor, mass), col, open_col)])
            rec[f"{col}_n_members_expected"] = len(member_states)
            rec[f"{col}_n_members_sensitive"] = sum(int(has_sens) for _, _, has_sens in member_states)
            rec[f"{col}_n_members_finite"] = sum(int(x is not None) for x, _, _ in member_states)
            rec[f"{col}_n_members_open"] = sum(int(is_open) for _, is_open, _ in member_states)
            topology_changed = any(
                (not has_sens) or (is_open != c_open) or ((x is not None) != (xc is not None))
                for x, is_open, has_sens in member_states)
            rec[f"{col}_topology_changed"] = topology_changed
            if xc is None:
                # central boundary open/undefined -> ribbon open on this side
                rec[f"{col}_band_lo"] = np.nan
                rec[f"{col}_band_hi"] = np.nan
                rec[f"{col}_open"] = True
                continue
            if topology_changed:
                rec[f"{col}_band_lo"] = np.nan
                rec[f"{col}_band_hi"] = np.nan
                rec[f"{col}_open"] = c_open or any(is_open for _, is_open, _ in member_states)
                continue
            any_open = c_open
            # --- scale: asymmetric envelope of central + scale curves ---
            xs = [xc]
            for name in axes.get("scale", []):
                x, op, _ = _boundary(curves[name], (flavor, mass), col, open_col)
                if op:
                    any_open = True
                if x is not None:
                    xs.append(x)
            scale_up = max(xs) - xc
            scale_dn = xc - min(xs)
            # --- pdf: std over replica members ---
            xp = []
            for name in axes.get("pdf", []):
                x, op, _ = _boundary(curves[name], (flavor, mass), col, open_col)
                if op:
                    any_open = True
                if x is not None:
                    xp.append(x)
            pdf_sigma = float(np.std(xp, ddof=1)) if len(xp) >= 2 else 0.0
            # --- mass: independent bottom/charm, per-quark max |dev|, in quadrature ---
            mdev = {"mb": 0.0, "mc": 0.0}
            for name in axes.get("mass", []):
                x, op, _ = _boundary(curves[name], (flavor, mass), col, open_col)
                if op:
                    any_open = True
                if x is None:
                    continue
                quark = "mb" if name.startswith("mb") else "mc"
                mdev[quark] = max(mdev[quark], abs(x - xc))
            mass_dev = float(np.hypot(mdev["mb"], mdev["mc"]))
            # --- alpha_s: half the |x| spread of the two companions, symmetric ---
            alphas_dev = 0.0
            if fold_alphas:
                xlo, op_lo, _ = _boundary(as_lo, (flavor, mass), col, open_col)
                xhi, op_hi, _ = _boundary(as_hi, (flavor, mass), col, open_col)
                if op_lo or op_hi:
                    any_open = True
                if xlo is not None and xhi is not None:
                    alphas_dev = 0.5 * abs(xhi - xlo)

            sigma_up = float(np.sqrt(scale_up**2 + pdf_sigma**2 + mass_dev**2 + alphas_dev**2))
            sigma_dn = float(np.sqrt(scale_dn**2 + pdf_sigma**2 + mass_dev**2 + alphas_dev**2))

            rec[f"{col}_band_lo"] = 10.0 ** (xc - sigma_dn)
            rec[f"{col}_band_hi"] = 10.0 ** (xc + sigma_up)
            rec[f"{col}_open"] = any_open
            rec[f"{col}_scale_up_dex"] = scale_up
            rec[f"{col}_scale_dn_dex"] = scale_dn
            rec[f"{col}_pdf_sigma_dex"] = pdf_sigma
            rec[f"{col}_mass_dev_dex"] = mass_dev
            if fold_alphas:
                rec[f"{col}_alphas_dex"] = alphas_dev
        if (flavor, float(mass)) in suppress:
            rec["u2_min_band_lo"] = np.nan
            rec["u2_min_band_hi"] = np.nan
            rec["band_suppressed"] = True
        out_rows.append(rec)
    out = pd.DataFrame(out_rows).sort_values(["flavor", "mass_GeV"])
    if "band_suppressed" in out.columns:
        out["band_suppressed"] = out["band_suppressed"].fillna(False).astype(bool)
    return out


def axis_members(registry: dict) -> dict[str, list[str]]:
    axes: dict[str, list[str]] = {axis: [] for axis in AXES}
    for v in registry["variations"]:
        axes[v["axis"]].append(v["name"])
    return axes


def combine_band(registry_path, alphas_lo=None, alphas_hi=None, suppress=frozenset()) -> pd.DataFrame:
    """Combine from a ``band_registry.json`` (one sensitivity CSV per variation)."""
    registry = json.loads(Path(registry_path).read_text())
    curves = {v["name"]: read_curve(v["sensitivity_csv"]) for v in registry["variations"]}
    alphas = None
    if alphas_lo and alphas_hi:
        alphas = (read_curve(alphas_lo), read_curve(alphas_hi))
    return combine_curves(curves, axis_members(registry), alphas=alphas, suppress=suppress)


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--registry", required=True, help="band_registry.json from the campaign")
    ap.add_argument("--out", default=None, help="output band CSV (default: <registry dir>/hnl_band.csv)")
    ap.add_argument("--alphas-lo", default=None,
                    help="as=0.117 companion sensitivity CSV; folded only with --alphas-hi")
    ap.add_argument("--alphas-hi", default=None, help="as=0.119 companion sensitivity CSV")
    ap.add_argument("--suppress", nargs="*", default=[], metavar="FLAVOR:MASS",
                    help="withhold the u2_min ribbon at these points (e.g. Ue:0.305)")
    args = ap.parse_args(argv)
    registry_path = Path(args.registry)
    if not registry_path.exists():
        print(f"registry not found: {registry_path}")
        return 1
    suppress = {(s.split(":")[0], float(s.split(":")[1])) for s in args.suppress}
    df = combine_band(registry_path, args.alphas_lo, args.alphas_hi, suppress)
    out = Path(args.out) if args.out else registry_path.parent / "hnl_band.csv"
    df.to_csv(out, index=False)
    axes = axis_members(json.loads(registry_path.read_text()))
    n_sens = int(df["has_sensitivity"].sum())
    note = ", alpha_s" if (args.alphas_lo and args.alphas_hi) else ""
    print(f"combined {sum(len(v) for v in axes.values())} variations "
          f"(scale={len(axes['scale'])}, pdf={len(axes['pdf'])}, mass={len(axes['mass'])}{note})")
    print(f"band written: {out}  ({n_sens}/{len(df)} mass points with sensitivity)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
