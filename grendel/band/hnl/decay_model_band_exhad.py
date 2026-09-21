"""Add the exHad decay model to the HNL decay-model band as a third member.

``decay_model_band`` varies the FairShip total width parametrically by
``delta(m)`` (5-20%) and re-scans the central run's acceptance Monte Carlo
with the lifetime and the hadronic/leptonic mix moved coherently; production
four-vectors and geometry are reused by construction. The exHad decay model
(arXiv:2609.16104) is an alternate model, not a rate factor: its width table
and its final states both differ from FairShip's.

This module adds exHad to that band the way the band itself was built -- on
the same production and geometry as its central, with its Monte-Carlo
settings (exact hit sample, 50 decay samples per hit, event chunk 1000) --
from two reruns of ``decay_model_band`` at the band's anchors:

* ``--control``: the FairShip templates. Where its ray-cast hit counts equal
  the band's, it is the band central recomputed by today's code; the
  agreement is recorded per anchor.
* ``--exhad``: the exHad templates.

The exHad member is the control-to-exHad shift in ``log10(U^2)`` applied to
the band's own central (the same dex rebasing the money plot uses), so
production and seed differences at anchors whose four-vectors changed since
the band was run do not enter the member.

Output: the band's rows and columns, with ``u2_*_dm_lo/hi`` now the ORDERED
span of the central, both width legs and the exHad member (NaN where any
member changes the topology); ``u2_*_width_lo/hi`` the original legs;
``u2_*_exhad`` the rebased member and ``u2_*_exhad_shift_dex`` its shift;
``u2_*_control``, ``u2_*_exhad50``, ``n_hits_control`` the rerun values.

    python -m grendel.band.hnl.decay_model_band_exhad --band decay_model_band.csv \\
        --control fairship50.csv --exhad exhad50.csv --out decay_model_band_exhad.csv
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ...io.atomic import sha256_file
from ..campaign_store import utc_now

KEYS = ["flavor", "mass_GeV"]
EDGES = ("u2_min", "u2_max")
REQUIRED_SETTINGS = {"decay_samples": 50, "event_chunk": 1000, "hit_estimator": "exact"}


def load(path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    for column in frame.columns:
        if column.startswith("has_sensitivity") or column.endswith("_open"):
            frame[column] = frame[column].astype(str).str.lower().isin(("true", "1"))
    return frame


def check_settings(frame: pd.DataFrame, label: str) -> None:
    for column, expected in REQUIRED_SETTINGS.items():
        values = set(frame[column].unique())
        if values != {expected}:
            raise RuntimeError(f"{label}: {column} = {sorted(values)}, the band requires {expected}")


def _edge(row, edge):
    """log10 of a finite, closed boundary of a sensitive row; None otherwise."""
    if row is None or not bool(row["has_sensitivity"]) or bool(row[f"{edge}_open"]):
        return None
    value = row[edge]
    return float(np.log10(value)) if np.isfinite(value) and value > 0 else None


def combine(band: pd.DataFrame, control: pd.DataFrame, exhad: pd.DataFrame):
    """The widened band and its bookkeeping."""
    control_rows = {(r.flavor, float(r.mass_GeV)): control.loc[(control.flavor == r.flavor) & (control.mass_GeV == r.mass_GeV)].iloc[0]
                    for r in control.itertuples(index=False)}
    exhad_rows = {(r.flavor, float(r.mass_GeV)): exhad.loc[(exhad.flavor == r.flavor) & (exhad.mass_GeV == r.mass_GeV)].iloc[0]
                  for r in exhad.itertuples(index=False)}
    rows, topology = [], {"control_without_central_boundary": [], "exhad_without_central_boundary": []}
    for _, b in band.sort_values(KEYS).iterrows():
        key = (b.flavor, float(b.mass_GeV))
        c, e = control_rows.get(key), exhad_rows.get(key)
        rec = b.to_dict()
        rec["n_hits_control"] = int(c["n_hits"]) if c is not None else -1
        rec["same_production_as_band"] = bool(c is not None and int(c["n_hits"]) == int(b["n_hits"]))
        for edge in EDGES:
            rec[f"{edge}_width_lo"] = b[f"{edge}_dm_lo"]
            rec[f"{edge}_width_hi"] = b[f"{edge}_dm_hi"]
            rec[f"{edge}_control"] = float(c[edge]) if c is not None else np.nan
            rec[f"{edge}_exhad50"] = float(e[edge]) if e is not None else np.nan
            xb, xc, xe = _edge(b, edge), _edge(c, edge), _edge(e, edge)
            rec[f"{edge}_control_ratio"] = 10.0 ** (xc - xb) if (xb is not None and xc is not None) else np.nan
            rec[f"{edge}_exhad_member_missing"] = bool(xb is not None and (xc is None or xe is None))
            if xb is not None and xc is None:
                topology["control_without_central_boundary"].append([b.flavor, key[1], edge])
            if xb is not None and xc is not None and xe is None:
                topology["exhad_without_central_boundary"].append([b.flavor, key[1], edge])
            if xb is not None and xc is not None and xe is not None:
                shift = xe - xc
                rec[f"{edge}_exhad_shift_dex"] = shift
                rec[f"{edge}_exhad"] = 10.0 ** (xb + shift)
            else:
                rec[f"{edge}_exhad_shift_dex"] = np.nan
                rec[f"{edge}_exhad"] = np.nan
            members = [b[edge], b[f"{edge}_dm_lo"], b[f"{edge}_dm_hi"], rec[f"{edge}_exhad"]]
            if xb is not None and all(np.isfinite(v) and v > 0 for v in members):
                rec[f"{edge}_dm_lo"] = float(min(members))
                rec[f"{edge}_dm_hi"] = float(max(members))
            else:
                # A member without this boundary is a topology change, not a ribbon.
                rec[f"{edge}_dm_lo"] = np.nan
                rec[f"{edge}_dm_hi"] = np.nan
        rows.append(rec)
    out = pd.DataFrame(rows)
    band_keys = set(zip(band.flavor, band.mass_GeV.astype(float)))
    topology["exhad_sensitive_where_central_is_not"] = [
        [r.flavor, float(r.mass_GeV), float(r.u2_min), float(r.u2_max)]
        for r in exhad.itertuples(index=False)
        if (r.flavor, float(r.mass_GeV)) not in band_keys and bool(r.has_sensitivity)]
    return out, topology


def summary(out: pd.DataFrame) -> dict:
    result = {"n_rows": int(len(out)), "n_same_production": int(out["same_production_as_band"].sum()),
              "different_production_anchors": [[r.flavor, float(r.mass_GeV), int(r.n_hits), int(r.n_hits_control)]
                                               for r in out[~out["same_production_as_band"]].itertuples(index=False)]}
    for edge in EDGES:
        ratio = out.loc[out["same_production_as_band"], f"{edge}_control_ratio"].dropna()
        shift = out[f"{edge}_exhad_shift_dex"].dropna()
        lo_old, hi_old = out[f"{edge}_width_lo"], out[f"{edge}_width_hi"]
        old_lo = np.fmin(np.fmin(lo_old, hi_old), out[edge])
        old_hi = np.fmax(np.fmax(lo_old, hi_old), out[edge])
        down = np.log10(old_lo / out[f"{edge}_dm_lo"])
        up = np.log10(out[f"{edge}_dm_hi"] / old_hi)
        finite = np.isfinite(down) & np.isfinite(up)
        result[edge] = {
            "control_vs_band_central_same_production": {
                "n": int(len(ratio)), "n_identical": int((ratio == 1.0).sum()),
                "max_abs_dex": float(np.abs(np.log10(ratio)).max()) if len(ratio) else None},
            "exhad_shift_dex": {"n": int(len(shift)), "median": float(shift.median()) if len(shift) else None,
                                "min": float(shift.min()) if len(shift) else None,
                                "max": float(shift.max()) if len(shift) else None},
            "band_with_interval": int(finite.sum()),
            "extended_down": int((down[finite] > 1e-12).sum()), "extended_up": int((up[finite] > 1e-12).sum()),
            "max_extension_down_dex": float(down[finite].max()) if finite.any() else None,
            "max_extension_up_dex": float(up[finite].max()) if finite.any() else None,
            "members_missing": int(out[f"{edge}_exhad_member_missing"].sum())}
        for flavor in ("Ue", "Umu", "Utau"):
            sel = out.loc[out.flavor == flavor, f"{edge}_exhad_shift_dex"].dropna()
            result[edge][f"exhad_shift_dex_{flavor}"] = (
                {"n": int(len(sel)), "median": float(sel.median()), "min": float(sel.min()), "max": float(sel.max())}
                if len(sel) else None)
    return result


def containment(curve: pd.DataFrame, band: pd.DataFrame, lo_col: str, hi_col: str, edge: str) -> dict:
    m = band.merge(curve, on=KEYS, suffixes=("_b", "_c"))
    m = m[m["has_sensitivity_b"].astype(bool) & m["has_sensitivity_c"].astype(bool)]
    value = m[f"{edge}_c"].to_numpy(float)
    lo, hi = m[lo_col].to_numpy(float), m[hi_col].to_numpy(float)
    ok = np.isfinite(value) & np.isfinite(lo) & np.isfinite(hi)
    value, lo, hi = value[ok], np.fmin(lo[ok], hi[ok]), np.fmax(lo[ok], hi[ok])
    inside = (value >= lo) & (value <= hi)
    excess = np.where(value > hi, np.log10(value / hi), np.log10(value / lo))
    excess = np.where(inside, 0.0, excess)
    return {"n": int(len(inside)), "inside": int(inside.sum()),
            "worst_breach_dex": float(np.abs(excess).max()) if len(excess) else None}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--band", required=True, help="decay_model_band.csv of the campaign")
    ap.add_argument("--control", required=True, help="decay_model_band rerun with the FairShip templates")
    ap.add_argument("--exhad", required=True, help="decay_model_band rerun with the exHad templates")
    ap.add_argument("--exhad-central", default=None, help="the exHad-decayed HNL curve, for the containment check")
    ap.add_argument("--out", required=True)
    ap.add_argument("--manifest-out", default=None, help="default: <out>.MANIFEST.json")
    ap.add_argument("--production", default="", help="production four-vector directory of the reruns (recorded)")
    ap.add_argument("--templates-control", default="", help="FairShip template directory (recorded)")
    ap.add_argument("--templates-exhad", default="", help="exHad template directory (recorded)")
    a = ap.parse_args(argv)
    band, control, exhad = load(a.band), load(a.control), load(a.exhad)
    for frame, label in ((band, "band"), (control, "control"), (exhad, "exhad")):
        check_settings(frame, label)
    if control.duplicated(KEYS).any() or exhad.duplicated(KEYS).any():
        raise RuntimeError("duplicate (flavor, mass) rows in a rerun")
    out, topology = combine(band, control, exhad)
    ordered = list(band.columns) + [c for c in out.columns if c not in band.columns]
    out = out.reindex(columns=ordered)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, index=False)
    result = summary(out)
    exhad_curve = load(a.exhad_central) if a.exhad_central and Path(a.exhad_central).exists() else None
    manifest = {
        "artifact": "GRENDEL HNL decay-model band, exHad member supplement",
        "schema_version": 1,
        "published_utc": utc_now(),
        "status": "supplement: decay_model_band.csv is unchanged; this file adds the exHad decay model "
                  "(arXiv:2609.16104) as a third member of the decay-model band, built like the band "
                  "itself on the same production and geometry with its Monte-Carlo settings",
        "construction": {
            "band_settings": REQUIRED_SETTINGS,
            "member": "log10(U^2) shift between the exHad and FairShip reruns on the campaign production, "
                      "applied to the band's own central at each anchor",
            "band": "u2_*_dm_lo/hi = ordered span of the central, the +/-delta legs (kept as "
                    "u2_*_width_lo/hi) and the exHad member; NaN where any member lacks the boundary",
            "production_reuse": "the decay-model band reuses the central run's four-vectors and geometry "
                                "by construction; so does this member"},
        "inputs": {Path(a.band).name: sha256_file(a.band),
                   "control_rerun": {"file": str(a.control), "sha256": sha256_file(a.control),
                                     "rows": int(len(control)), "templates": a.templates_control},
                   "exhad_rerun": {"file": str(a.exhad), "sha256": sha256_file(a.exhad),
                                   "rows": int(len(exhad)), "templates": a.templates_exhad},
                   "production": a.production},
        "summary": result,
        "topology": topology,
        "exhad_central_curve_containment": (
            {"curve": Path(a.exhad_central).name, "sha256": sha256_file(a.exhad_central),
             **{edge: {"band": containment(exhad_curve, band, f"{edge}_dm_lo", f"{edge}_dm_hi", edge),
                       "supplemented_band": containment(exhad_curve, out, f"{edge}_dm_lo", f"{edge}_dm_hi", edge)}
                for edge in EDGES}} if exhad_curve is not None else None),
        "outputs": {Path(a.out).name: {"sha256": sha256_file(a.out), "rows": int(len(out))}},
    }
    manifest_out = Path(a.manifest_out) if a.manifest_out else Path(a.out).with_suffix(".MANIFEST.json")
    manifest_out.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {a.out} ({len(out)} rows) and {manifest_out}")
    print(json.dumps({"summary": result, "topology": topology,
                      "containment": manifest["exhad_central_curve_containment"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
