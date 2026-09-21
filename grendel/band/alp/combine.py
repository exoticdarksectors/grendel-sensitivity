"""Combine the BC10 campaign's variation curves into the display envelope.

Named alternatives (scale, bottom mass, gluon surrogate, C_bs) use their
extrema together with the central; the NNPDF ensemble uses its 16th and
84th percentiles. The display envelope is the outermost boundary of those
one-source intervals in ``log10(1/f)``; same-physics repeats are numerical
controls, reported but excluded. No sources are combined in quadrature and
this is not a confidence interval (see :mod:`grendel.band.envelope`).

The campaign measures one-source shifts relative to its own nominal
sample. :func:`combine_band` optionally rebases the pointwise halo
components onto a canonical, higher-statistics central contour; the
numerical repeats keep their directly simulated values.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..envelope import axis_interval, log10_boundary, pdf_percentiles, repeat_statistics

BOUNDARIES = ("invf_min", "invf_max")
HALO_AXES = (("scale", "scale"), ("mb", "mb"), ("decay_gg", "gg"), ("cbs", "cbs"))
EXPECTED_COUNTS = {"central": 1, "scale": 6, "pdf": 100, "mb": 2,
                   "decay_gg": 3, "cbs": 2, "numerical_control": 2}


def _log_boundary(raw, mass, variation, boundary):
    """(x, is_open, sensitive) for one variation at one mass in the raw table."""
    selected = raw[(raw["mass_GeV"] == mass) & (raw["variation"] == variation)]
    if len(selected) != 1:
        return None, False, False
    row = selected.iloc[0]
    sensitive = bool(row["has_sensitivity"])
    opened = bool(row.get(f"{boundary}_open", False))
    if not sensitive:
        return None, opened, sensitive
    return log10_boundary(row.get(boundary, np.nan)), opened, sensitive


def combine_campaign_band(raw: pd.DataFrame) -> pd.DataFrame:
    """The pointwise one-source-at-a-time envelope in campaign coordinates."""
    axes = {axis: list(raw.loc[raw["axis"] == axis, "variation"].unique())
            for axis in EXPECTED_COUNTS}
    counts = {axis: len(names) for axis, names in axes.items()}
    if counts != EXPECTED_COUNTS:
        raise ValueError(f"incomplete contour ensemble: {counts} != {EXPECTED_COUNTS}")
    halo_variations = [name for axis, _ in HALO_AXES for name in axes[axis]] + axes["pdf"]
    halo_variations = [name for axis in ("scale", "pdf", "mb", "decay_gg", "cbs") for name in axes[axis]]
    rows = []
    for mass in sorted(raw["mass_GeV"].unique()):
        central = raw[(raw["mass_GeV"] == mass) & (raw["variation"] == "central")].iloc[0]
        central_sensitive = bool(central["has_sensitivity"])
        numerical_rows = raw[(raw["mass_GeV"] == mass)
                             & raw["variation"].isin(axes["numerical_control"])].set_index("variation")
        numerical_sensitive = {name: bool(numerical_rows.loc[name, "has_sensitivity"])
                               for name in axes["numerical_control"]}
        numerical_topology_differences = {
            name for name, sensitive in numerical_sensitive.items() if sensitive != central_sensitive}
        halo_rows = raw[(raw["mass_GeV"] == mass)
                        & raw["variation"].isin(halo_variations)].set_index("variation")
        halo_sensitive = {name: bool(halo_rows.loc[name, "has_sensitivity"]) for name in halo_variations}
        halo_restores_sensitivity = {name for name, s in halo_sensitive.items() if s and not central_sensitive}
        halo_removes_sensitivity = {name for name, s in halo_sensitive.items() if central_sensitive and not s}
        halo_topology_differences = halo_restores_sensitivity | halo_removes_sensitivity
        any_halo_variation_sensitive = bool(
            raw.loc[(raw["mass_GeV"] == mass) & ~raw["axis"].isin(["numerical_control"]),
                    "has_sensitivity"].any())
        record = {
            "mass_GeV": mass,
            "has_sensitivity": central_sensitive,
            "any_variation_sensitive": any_halo_variation_sensitive,
            "any_halo_variation_sensitive": any_halo_variation_sensitive,
            "numerical_control_n_sensitive": sum(numerical_sensitive.values()),
            "numerical_control_any_sensitive": any(numerical_sensitive.values()),
            "numerical_control_all_sensitive": all(numerical_sensitive.values()),
            "numerical_control_sensitive_variations": ";".join(
                name for name in axes["numerical_control"] if numerical_sensitive[name]),
            "numerical_control_included_in_halo": False,
            "halo_restores_sensitivity": bool(halo_restores_sensitivity),
            "halo_restores_sensitivity_variations": ";".join(
                name for name in halo_variations if name in halo_restores_sensitivity),
            "halo_removes_sensitivity": bool(halo_removes_sensitivity),
            "halo_removes_sensitivity_variations": ";".join(
                name for name in halo_variations if name in halo_removes_sensitivity),
            "envelope_definition": "single_source_variation_envelope",
        }
        for boundary in BOUNDARIES:
            xc, central_open, _ = _log_boundary(raw, mass, "central", boundary)
            record[f"{boundary}_central"] = float(central.get(boundary, np.nan))
            record[f"{boundary}_open"] = central_open
            for name in halo_variations:
                _, opened, sensitive = _log_boundary(raw, mass, name, boundary)
                if central_sensitive and sensitive and central_open != opened:
                    halo_topology_differences.add(name)
            repeat_samples = []
            for name in axes["numerical_control"]:
                value, opened, sensitive = _log_boundary(raw, mass, name, boundary)
                repeat_samples.append((name, value, opened, sensitive))
                record[f"{boundary}_{name}"] = 10.0 ** value if value is not None else np.nan
                if central_sensitive and sensitive and central_open != opened:
                    numerical_topology_differences.add(name)
            if xc is None:
                record.update({
                    f"{boundary}_envelope_lo": np.nan, f"{boundary}_envelope_hi": np.nan,
                    f"{boundary}_any_variation_open": central_open,
                    f"{boundary}_variation_missing": False,
                    f"{boundary}_envelope_lo_source": "", f"{boundary}_envelope_hi_source": "",
                    f"{boundary}_physical_envelope_max_abs_dex": np.nan,
                    f"{boundary}_repeat_median_abs_dex": np.nan, f"{boundary}_repeat_max_abs_dex": np.nan,
                    f"{boundary}_repeat_median_abs_fraction": np.nan,
                    f"{boundary}_repeat_max_abs_fraction": np.nan,
                    f"{boundary}_repeat_missing": False, f"{boundary}_repeat_not_subdominant": False,
                })
                continue
            missing = False
            any_open = central_open

            def values(names):
                nonlocal missing, any_open
                result = []
                for name in names:
                    value, opened, _ = _log_boundary(raw, mass, name, boundary)
                    any_open = any_open or opened
                    if value is None:
                        missing = True
                    else:
                        result.append((name, value))
                return result

            source_intervals = {}
            for axis, label in HALO_AXES:
                interval = axis_interval(values(axes[axis]), central=("central", xc))
                source_intervals[label] = (interval.lo, interval.hi)
                record.update({
                    f"{boundary}_{label}_lo": 10.0 ** interval.lo,
                    f"{boundary}_{label}_hi": 10.0 ** interval.hi,
                    f"{boundary}_{label}_lo_variation": interval.lo_name,
                    f"{boundary}_{label}_hi_variation": interval.hi_name,
                })
            pdf = pdf_percentiles(value for _, value in values(axes["pdf"]))
            pdf_p16, pdf_p84 = (pdf.p16, pdf.p84) if pdf.n else (xc, xc)
            source_intervals["pdf"] = (min(xc, float(pdf_p16)), max(xc, float(pdf_p84)))
            record.update({
                f"{boundary}_pdf_p16": 10.0 ** float(pdf_p16),
                f"{boundary}_pdf_p84": 10.0 ** float(pdf_p84),
                f"{boundary}_pdf_log10_std": pdf.std,
                f"{boundary}_pdf_n_finite": pdf.n,
            })
            lo_source, lo = min(((s, iv[0]) for s, iv in source_intervals.items()), key=lambda item: item[1])
            hi_source, hi = max(((s, iv[1]) for s, iv in source_intervals.items()), key=lambda item: item[1])
            record.update({
                f"{boundary}_envelope_lo": 10.0 ** lo,
                f"{boundary}_envelope_hi": 10.0 ** hi,
                f"{boundary}_any_variation_open": any_open,
                f"{boundary}_variation_missing": missing,
                f"{boundary}_envelope_lo_source": lo_source,
                f"{boundary}_envelope_hi_source": hi_source,
            })
            repeat_values = [value for _, value, _, _ in repeat_samples if value is not None]
            repeat_missing = any(value is None for _, value, _, _ in repeat_samples)
            repeat = repeat_statistics(xc, repeat_values)
            physical_max_abs_dex = max(xc - lo, hi - xc)
            record.update({
                f"{boundary}_physical_envelope_max_abs_dex": physical_max_abs_dex,
                f"{boundary}_repeat_median_abs_dex": repeat.median_abs_dex,
                f"{boundary}_repeat_max_abs_dex": repeat.max_abs_dex,
                f"{boundary}_repeat_median_abs_fraction": repeat.median_abs_fraction,
                f"{boundary}_repeat_max_abs_fraction": repeat.max_abs_fraction,
                f"{boundary}_repeat_missing": repeat_missing,
                f"{boundary}_repeat_not_subdominant": bool(
                    repeat_missing
                    or (np.isfinite(repeat.max_abs_dex) and repeat.max_abs_dex > 0.0
                        and repeat.max_abs_dex >= physical_max_abs_dex)),
            })
        record["halo_topology_differs"] = bool(halo_topology_differences)
        record["halo_topology_difference_variations"] = ";".join(
            name for name in halo_variations if name in halo_topology_differences)
        record["numerical_control_topology_differs"] = bool(numerical_topology_differences)
        record["numerical_control_topology_difference_variations"] = ";".join(
            name for name in axes["numerical_control"] if name in numerical_topology_differences)
        rows.append(record)
    return pd.DataFrame(rows).sort_values("mass_GeV")


def combine_band(raw: pd.DataFrame, central_curve=None) -> pd.DataFrame:
    """Build the campaign band and optionally rebase it to a canonical contour.

    Only the pointwise halo components are rebased (multiplied by
    canonical/campaign at each mass); the numerical repeats retain their
    independently simulated absolute values.
    """
    band = combine_campaign_band(raw)
    if central_curve is None:
        return band
    reference = (pd.read_csv(central_curve) if isinstance(central_curve, (str, Path))
                 else central_curve.copy())
    if reference["mass_GeV"].duplicated().any():
        raise ValueError("canonical central curve contains duplicate masses")
    reference = reference.set_index("mass_GeV")
    missing = sorted(set(band["mass_GeV"]) - set(reference.index.astype(float)))
    if missing:
        raise ValueError(f"canonical central curve is missing {len(missing)} campaign masses; "
                         f"first: {missing[0]}")
    halo_value_columns = {
        boundary: [*(f"{boundary}_{label}_{edge}" for _, label in HALO_AXES for edge in ("lo", "hi")),
                   f"{boundary}_pdf_p16", f"{boundary}_pdf_p84",
                   f"{boundary}_envelope_lo", f"{boundary}_envelope_hi"]
        for boundary in BOUNDARIES
    }
    for index, campaign in band.iterrows():
        mass = float(campaign["mass_GeV"])
        canonical = reference.loc[mass]
        campaign_sensitive = bool(campaign["has_sensitivity"])
        canonical_sensitive = bool(canonical.get("has_sensitivity", False))
        band.at[index, "campaign_has_sensitivity"] = campaign_sensitive
        band.at[index, "has_sensitivity"] = canonical_sensitive
        band.at[index, "any_variation_sensitive"] = bool(campaign["any_variation_sensitive"] or canonical_sensitive)
        band.at[index, "any_halo_variation_sensitive"] = bool(
            campaign["any_halo_variation_sensitive"] or canonical_sensitive)
        topology_compatible = campaign_sensitive == canonical_sensitive
        for boundary in BOUNDARIES:
            campaign_value = float(campaign[f"{boundary}_central"])
            canonical_value = float(canonical.get(boundary, np.nan))
            campaign_open = bool(campaign[f"{boundary}_open"])
            canonical_open = bool(canonical.get(f"{boundary}_open", False))
            band.at[index, f"{boundary}_campaign_central"] = campaign_value
            band.at[index, f"{boundary}_central"] = canonical_value
            band.at[index, f"{boundary}_open"] = canonical_open
            compatible = (campaign_sensitive and canonical_sensitive
                          and np.isfinite(campaign_value) and campaign_value > 0.0
                          and np.isfinite(canonical_value) and canonical_value > 0.0
                          and campaign_open == canonical_open and not campaign_open)
            topology_compatible = topology_compatible and (compatible or not campaign_sensitive)
            if compatible:
                factor = canonical_value / campaign_value
                for column in halo_value_columns[boundary]:
                    value = band.at[index, column]
                    if np.isfinite(value):
                        band.at[index, column] = float(value) * factor
            else:
                for column in halo_value_columns[boundary]:
                    band.at[index, column] = np.nan
                band.at[index, f"{boundary}_variation_missing"] = True
        band.at[index, "canonical_rebase_topology_compatible"] = bool(topology_compatible)
    band["envelope_reference"] = "canonical_high_statistics_central"
    return band
