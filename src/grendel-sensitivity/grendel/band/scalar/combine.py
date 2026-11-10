"""Combine the BC4 campaign's variation curves into the display envelope.

Scale and bottom-mass sources use their standard extrema, PDFs use the
16th/84th percentiles of the 100 replica boundaries, and the decay source
is the interval spanned by the central and every ``decay_model``-axis
alternate (the LO-ChPT/spectator run, and the exHad supplement when
present). The display envelope is the outermost endpoint among those
one-source-at-a-time intervals; every interval is rebased in dex onto the
canonical central curve. Same-physics repeats are numerical controls,
reported but excluded. It is neither a quadrature combination nor a
confidence band (see :mod:`grendel.band.envelope`).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..envelope import axis_interval, log10_boundary, pdf_percentiles, repeat_statistics

BOUNDARIES = (("u2_min", "u2_min_open"), ("u2_max", "u2_max_open"))
DECAY_VARIATION = "decay_chpt_spectator"
PHYSICS_AXES = ("scale", "pdf", "mass", "decay_model")


def _boundary(row, boundary, open_col):
    """(x, is_open) of a sensitive row's boundary; (None, False) when not sensitive."""
    if row is None or not bool(row.get("has_sensitivity", False)):
        return None, False
    is_open = bool(row.get(open_col, False))
    return log10_boundary(row.get(boundary, np.nan)), is_open


def combine_band(raw, central_curve) -> pd.DataFrame:
    """The envelope around ``central_curve`` from the raw variation table."""
    raw = pd.read_csv(raw) if isinstance(raw, (str, Path)) else raw.copy()
    reference = (pd.read_csv(central_curve) if isinstance(central_curve, (str, Path))
                 else central_curve.copy())
    axes = {axis: sorted(raw.loc[raw["axis"] == axis, "variation"].unique())
            for axis in (*PHYSICS_AXES, "numerical_control")}
    physical_variations = [name for axis in PHYSICS_AXES for name in axes[axis]]
    rows = []
    for _, ref in reference.sort_values("mass_GeV").iterrows():
        mass = float(ref["mass_GeV"])
        group = raw[np.isclose(raw["mass_GeV"], mass, rtol=0, atol=5e-10)]
        by_name = {row["variation"]: row for _, row in group.iterrows()}
        campaign_central = by_name.get("central")
        physics_group = group[group["axis"] != "numerical_control"]
        campaign_sensitive = bool(campaign_central is not None
                                  and campaign_central.get("has_sensitivity", False))
        repeat_sensitive = {
            name: bool(by_name.get(name) is not None and by_name[name].get("has_sensitivity", False))
            for name in axes["numerical_control"]}
        topology_differs_from_campaign = {
            name for name, sensitive in repeat_sensitive.items() if sensitive != campaign_sensitive}
        canonical_sensitive = bool(ref.get("has_sensitivity", False))
        topology_differs_from_canonical = {
            name for name, sensitive in repeat_sensitive.items() if sensitive != canonical_sensitive}
        physical_sensitive = {
            name: bool(by_name.get(name) is not None and by_name[name].get("has_sensitivity", False))
            for name in physical_variations}
        physical_topology_differs_from_campaign = {
            name for name, sensitive in physical_sensitive.items() if sensitive != campaign_sensitive}
        physical_topology_differs_from_canonical = {
            name for name, sensitive in physical_sensitive.items() if sensitive != canonical_sensitive}
        rec = {
            "mass_GeV": mass,
            "has_sensitivity": canonical_sensitive,
            "any_variation_sensitive": bool(
                len(physics_group)
                and physics_group["has_sensitivity"].fillna(False).astype(bool).any()),
            "envelope_definition": "single_source_variation_envelope",
            "campaign_has_sensitivity": campaign_sensitive,
            "numerical_control_n_sensitive": sum(repeat_sensitive.values()),
            "numerical_control_any_sensitive": any(repeat_sensitive.values()),
            "numerical_control_all_sensitive": bool(repeat_sensitive) and all(repeat_sensitive.values()),
            "numerical_control_sensitive_variations": ";".join(
                name for name in axes["numerical_control"] if repeat_sensitive[name]),
            "numerical_control_included_in_envelope": False,
            "physical_variation_n_sensitive": sum(physical_sensitive.values()),
            "physical_variation_any_sensitive": any(physical_sensitive.values()),
            "physical_variation_sensitive_variations": ";".join(
                name for name in physical_variations if physical_sensitive[name]),
        }
        for boundary, open_col in BOUNDARIES:
            ref_value = ref.get(boundary, np.nan)
            ref_open = bool(ref.get(open_col, False))
            rec[f"{boundary}_central"] = ref_value
            rec[f"{boundary}_open"] = ref_open
            xc, campaign_open = _boundary(campaign_central, boundary, open_col)
            rec[f"{boundary}_campaign_central"] = 10.0**xc if xc is not None else np.nan
            rec[f"{boundary}_campaign_central_open"] = campaign_open
            repeat_samples = []
            for name in axes["numerical_control"]:
                value, is_open = _boundary(by_name.get(name), boundary, open_col)
                repeat_samples.append((name, value, is_open))
                rec[f"{boundary}_{name}"] = 10.0**value if value is not None else np.nan
                if campaign_sensitive and repeat_sensitive[name] and campaign_open != is_open:
                    topology_differs_from_campaign.add(name)
                if canonical_sensitive and repeat_sensitive[name] and ref_open != is_open:
                    topology_differs_from_canonical.add(name)
            for name in physical_variations:
                _, is_open = _boundary(by_name.get(name), boundary, open_col)
                if campaign_sensitive and physical_sensitive[name] and campaign_open != is_open:
                    physical_topology_differs_from_campaign.add(name)
                if canonical_sensitive and physical_sensitive[name] and ref_open != is_open:
                    physical_topology_differs_from_canonical.add(name)
            if (not rec["has_sensitivity"] or ref_open or not np.isfinite(ref_value)
                    or ref_value <= 0 or xc is None or campaign_open):
                rec[f"{boundary}_envelope_lo"] = np.nan
                rec[f"{boundary}_envelope_hi"] = np.nan
                rec[f"{boundary}_envelope_open"] = ref_open or campaign_open
                continue
            xref = float(np.log10(ref_value))
            any_open = ref_open or campaign_open
            missing = False

            def collect_axis(axis):
                nonlocal any_open, missing
                values = []
                for name in axes[axis]:
                    value, is_open = _boundary(by_name.get(name), boundary, open_col)
                    any_open |= is_open
                    missing |= value is None and not is_open
                    if value is not None and not is_open:
                        values.append((name, value))
                return values

            def rebased(value):
                return 10.0**(xref + value - xc)

            scale = axis_interval(collect_axis("scale"), central=("central", xc))
            pdf_named = collect_axis("pdf")
            pdf = pdf_percentiles(value for _, value in pdf_named)
            if pdf.n == 0:
                missing = True
            mass_iv = axis_interval(collect_axis("mass"), central=("central", xc))
            # Every decay-model alternate spans the interval together with the
            # central; the LO-ChPT/spectator run keeps its legacy columns.
            decay_named = collect_axis("decay_model")
            decay_by_name = {name: value for name, value in decay_named}
            alt = decay_by_name.get(DECAY_VARIATION)
            decay = axis_interval(decay_named, central=("central", xc))
            repeat_values, repeat_open, repeat_missing = [], False, False
            for name, value, is_open in repeat_samples:
                repeat_open |= is_open
                repeat_missing |= value is None and not is_open
                if value is not None and not is_open:
                    repeat_values.append(value)
            repeat = repeat_statistics(xc, repeat_values)
            physical_magnitudes = {
                "scale": scale.magnitude(xc),
                "pdf": max(abs(pdf.p16 - xc), abs(pdf.p84 - xc)),
                "bottom_mass": mass_iv.magnitude(xc),
                "decay_model": decay.magnitude(xc),
            }

            def numerical_ratio(component):
                magnitude = physical_magnitudes[component]
                if not np.isfinite(repeat.max_abs_dex) or magnitude <= 0:
                    return np.nan
                return repeat.max_abs_dex / magnitude

            display_magnitude = max(physical_magnitudes.values())
            numerical_subdominant = (
                bool(repeat.max_abs_dex < display_magnitude)
                if np.isfinite(repeat.max_abs_dex) and display_magnitude > 0 else False)
            source_endpoints = [
                (scale.lo, f"scale:{scale.lo_name}"), (scale.hi, f"scale:{scale.hi_name}"),
                (mass_iv.lo, f"bottom_mass:{mass_iv.lo_name}"), (mass_iv.hi, f"bottom_mass:{mass_iv.hi_name}"),
                (decay.lo, f"decay_model:{decay.lo_name}"), (decay.hi, f"decay_model:{decay.hi_name}"),
            ]
            if np.isfinite(pdf.p16) and np.isfinite(pdf.p84):
                source_endpoints.extend([(pdf.p16, "pdf:p16"), (pdf.p84, "pdf:p84")])
            total_lo = min(source_endpoints, key=lambda item: item[0])
            total_hi = max(source_endpoints, key=lambda item: item[0])
            rec.update({
                f"{boundary}_scale_envelope_lo": rebased(scale.lo),
                f"{boundary}_scale_envelope_hi": rebased(scale.hi),
                f"{boundary}_scale_envelope_lo_source": scale.lo_name,
                f"{boundary}_scale_envelope_hi_source": scale.hi_name,
                f"{boundary}_scale_envelope_lo_dex": scale.lo - xc,
                f"{boundary}_scale_envelope_hi_dex": scale.hi - xc,
                f"{boundary}_pdf_p16": rebased(pdf.p16) if np.isfinite(pdf.p16) else np.nan,
                f"{boundary}_pdf_p84": rebased(pdf.p84) if np.isfinite(pdf.p84) else np.nan,
                f"{boundary}_pdf_p16_shift_dex": pdf.p16 - xc if np.isfinite(pdf.p16) else np.nan,
                f"{boundary}_pdf_p84_shift_dex": pdf.p84 - xc if np.isfinite(pdf.p84) else np.nan,
                f"{boundary}_pdf_std_dex": pdf.std,
                f"{boundary}_pdf_n_finite": pdf.n,
                f"{boundary}_bottom_mass_envelope_lo": rebased(mass_iv.lo),
                f"{boundary}_bottom_mass_envelope_hi": rebased(mass_iv.hi),
                f"{boundary}_bottom_mass_envelope_lo_source": mass_iv.lo_name,
                f"{boundary}_bottom_mass_envelope_hi_source": mass_iv.hi_name,
                f"{boundary}_bottom_mass_envelope_lo_dex": mass_iv.lo - xc,
                f"{boundary}_bottom_mass_envelope_hi_dex": mass_iv.hi - xc,
                f"{boundary}_decay_model_alt": rebased(alt) if alt is not None else np.nan,
                f"{boundary}_decay_model_shift_dex": alt - xc if alt is not None else np.nan,
                **{f"{boundary}_decay_model_{name}": (
                    rebased(decay_by_name[name]) if name in decay_by_name else np.nan)
                   for name in axes["decay_model"]},
                **{f"{boundary}_decay_model_{name}_shift_dex": (
                    decay_by_name[name] - xc if name in decay_by_name else np.nan)
                   for name in axes["decay_model"]},
                f"{boundary}_decay_model_envelope_lo": rebased(decay.lo),
                f"{boundary}_decay_model_envelope_hi": rebased(decay.hi),
                f"{boundary}_decay_model_envelope_lo_source": decay.lo_name,
                f"{boundary}_decay_model_envelope_hi_source": decay.hi_name,
                f"{boundary}_numerical_repeat_median_abs_dex": repeat.median_abs_dex,
                f"{boundary}_numerical_repeat_max_abs_dex": repeat.max_abs_dex,
                f"{boundary}_numerical_repeat_max_fractional": repeat.max_abs_fraction,
                f"{boundary}_numerical_repeat_n_finite": repeat.n,
                f"{boundary}_numerical_repeat_open": repeat_open,
                f"{boundary}_numerical_repeat_missing": repeat_missing,
                f"{boundary}_numerical_repeat_to_scale_ratio": numerical_ratio("scale"),
                f"{boundary}_numerical_repeat_to_pdf_ratio": numerical_ratio("pdf"),
                f"{boundary}_numerical_repeat_to_bottom_mass_ratio": numerical_ratio("bottom_mass"),
                f"{boundary}_numerical_repeat_to_decay_model_ratio": numerical_ratio("decay_model"),
                f"{boundary}_numerical_repeat_subdominant": numerical_subdominant,
                f"{boundary}_numerical_repeat_not_subdominant": not numerical_subdominant,
                f"{boundary}_envelope_lo": rebased(total_lo[0]),
                f"{boundary}_envelope_hi": rebased(total_hi[0]),
                f"{boundary}_envelope_lo_source": total_lo[1],
                f"{boundary}_envelope_hi_source": total_hi[1],
                f"{boundary}_envelope_open": any_open,
                f"{boundary}_variation_missing": missing,
            })
        rec["numerical_control_topology_differs_from_campaign"] = bool(topology_differs_from_campaign)
        rec["numerical_control_topology_difference_variations_from_campaign"] = ";".join(
            name for name in axes["numerical_control"] if name in topology_differs_from_campaign)
        rec["numerical_control_topology_differs_from_canonical"] = bool(topology_differs_from_canonical)
        rec["numerical_control_topology_difference_variations_from_canonical"] = ";".join(
            name for name in axes["numerical_control"] if name in topology_differs_from_canonical)
        rec["physical_variation_topology_differs_from_campaign"] = bool(physical_topology_differs_from_campaign)
        rec["physical_variation_topology_difference_variations_from_campaign"] = ";".join(
            name for name in physical_variations if name in physical_topology_differs_from_campaign)
        rec["physical_variation_topology_differs_from_canonical"] = bool(physical_topology_differs_from_canonical)
        rec["physical_variation_topology_difference_variations_from_canonical"] = ";".join(
            name for name in physical_variations if name in physical_topology_differs_from_canonical)
        rows.append(rec)
    return pd.DataFrame(rows)
