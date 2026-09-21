"""BC4 / BC10 panels: existing exclusions overlaid with the GRENDEL island.

Gray filled existing-exclusion regions, the red GRENDEL island band, a
run-label box, log-log axes.

* BC4 (Higgs-portal dark scalar, sin^2 theta vs m_S): excluded region from
  the PBC report BC5 figure (mixing-only limits, squared from sin theta;
  ``curves.tools.vector_excluded_bc4``); SHiP and CODEX-b projections.
* BC10 (fermiophilic ALP, 1/f vs m_a, BNT convention c_f = 1): LHCb
  B -> K mu mu recast of GKOZ arXiv:2310.03524 (``curves.tools.vector_excluded_bc10``),
  the boundary chains clipped at the source frame top; past beam dumps and
  the FIPs 2022 bounds; SHiP and DarkQuest projections.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np

from .common import (GRENDEL_DARK_RED, GRENDEL_RED, SIGNAL_THRESHOLD, load_polygon_segments,
                     parse_csv_bool, parse_csv_float)
from .data import BC4_DATA, BC10_DATA

LEGEND_TITLE = "Curve data sources"

# Competitor projections: dedicated LLP experiments only. Sources and
# physics-validity notes live with the digitisation tools.
BENCHMARKS = {
    "bc4": {
        "prefix": "u2",
        "excluded": BC4_DATA / "vector" / "currently_excluded.dat",
        "excluded_closed": True,
        "excluded_label": "Currently excluded\narXiv:2505.00947",
        "extra_excluded": [],
        "projections": [
            {"path": BC4_DATA / "projections" / "ship_bc4.dat",
             "label": "SHiP\narXiv:2310.17726", "color": "#1f77b4", "ls": (0, (5, 2.2)), "x_clip": None},
            {"path": BC4_DATA / "projections" / "codexb_300fb.dat",
             "label": "CODEX-b\narXiv:1911.00481", "color": "#2ca02c", "ls": (0, (5, 1.6, 1, 1.6)),
             "x_clip": None},
            # A BC5 curve on a BC4 axis: FASER2 cannot produce Higgs bosons, so
            # BC5's BR(h->SS) is irrelevant to it -- but BC5's quartic also
            # opens B -> K S S, which a forward B-driven detector does see.
            # Measured on SHiP the BC5 lower branch is 10-60x deeper below
            # ~4 GeV, so this curve overstates FASER2's BC4 reach by an order
            # of magnitude at small mixing. Kept as a labelled diagnostic
            # only; never in the paper profile.
            {"path": BC4_DATA / "projections" / "faser2_pbc_bc5.dat",
             "label": "FASER2 (BC5)\narXiv:2505.00947", "color": "#9467bd", "ls": (0, (1.5, 1.8)),
             "x_clip": None},
        ],
        "xlabel": r"$m_S$ (GeV)",
        "ylabel": r"$\sin^2\theta$",
        # The closed island begins at the dimuon threshold 2 m_mu = 0.211 GeV
        # (first closed grid point 0.22 GeV); grendel_m_min drops the open-
        # topped e e -only rows below it. The view starts a little above
        # threshold to trim the very narrow tip where the upper edge overshoots.
        "x_range": (0.25, 5.0),
        "x_ticks": (0.3, 0.5, 1.0, 2.0, 5.0),
        "grendel_m_min": 0.211,
        # y-floor below the SHiP lower tip at 7.9e-13 so the whole SHiP curve
        # is shown; GRENDEL bottoms at 6.5e-12, well inside the frame.
        "y_range": (1.0e-13, 1.0e-4),
        "scan_ceiling": 1.0e-2,
        "legend_loc": "lower left",
        "run_label_tag": "PBC BC4 (dark scalar)",
        "figure": "bc4_grendel_talk_zoom",
        "x_clip": None,
    },
    "bc10": {
        "prefix": "invf",
        # Suppress islands narrower than 50 MeV: the two-point 1.40-1.41 GeV
        # crossing (peak_N 3.19/3.16, 6% over threshold, absent from the
        # N >= 10 contour) drew as a detached speck between the 1.25 and
        # 1.47 GeV lobes. A presentation choice, not a correction: the pocket
        # is real (an independent six-million-event control at 1.40 GeV
        # reproduces it). Anything removed is announced on stderr.
        "min_island_GeV": 0.05,
        "excluded": BC10_DATA / "vector" / "lhcb_bkmumu_excluded.dat",
        "excluded_closed": False,
        "excluded_label": "LHCb $B\\!\\to\\! K\\mu\\mu$ excluded\nGKOZ arXiv:2310.03524",
        "extra_excluded": [
            {"path": BC10_DATA / "vector" / "beamdumps_past_alp2.dat",
             "label": "Past beam dumps\narXiv:2501.04525", "mode": "closed"},
            {"path": BC10_DATA / "vector" / "na62_kpia_excluded.dat",
             "label": "NA62 $K^+\\!\\to\\!\\pi^+ + \\mathrm{inv}$\nvia ALPINIST", "mode": "band"},
            {"path": BC10_DATA / "vector" / "fips2022" / "kpix_upper.dat",
             "label": "FIPs 2022 existing bounds\n(pre-GKOZ phenomenology)", "mode": "closed"},
            {"path": BC10_DATA / "vector" / "fips2022" / "e949_kpi_inv.dat", "label": None, "mode": "closed"},
            {"path": BC10_DATA / "vector" / "fips2022" / "kl_pi0ll.dat", "label": None, "mode": "closed"},
            {"path": BC10_DATA / "vector" / "fips2022" / "b_to_k_inv.dat", "label": None, "mode": "closed"},
            {"path": BC10_DATA / "vector" / "fips2022" / "babar.dat", "label": None, "mode": "closed"},
            {"path": BC10_DATA / "vector" / "fips2022" / "bs_mumu.dat", "label": None, "mode": "closed"},
            {"path": BC10_DATA / "vector" / "fips2022" / "bbn.dat", "label": None, "mode": "closed"},
        ],
        "projections": [
            {"path": BC10_DATA / "projections" / "ship_alp2.dat",
             "label": "SHiP\narXiv:2501.04525", "color": "#1f77b4", "ls": (0, (5, 2.2)),
             "x_clip": (0.6, 2.0)},   # published figure frame; the .dat keeps the clipped path too
            {"path": BC10_DATA / "projections" / "darkquest_alp2.dat",
             "label": "DarkQuest phase I\narXiv:2501.04525", "color": "#e377c2",
             "ls": (0, (5, 1.6, 1, 1.6)), "x_clip": (0.6, 2.0)},
        ],
        "xlabel": r"$m_a$ (GeV)",
        # The reference curves store 1/f in the BNT convention:
        # g_Y = v_h/f_BNT = 2 v_h/f_GKOZ, so the plotted scale factor is v_h.
        "ylabel": r"$g_Y = v_h/f$",
        "y_scale": 246.0,
        "secondary_ylabel": r"$1/f$ (GeV$^{-1}$)",
        "x_range": (0.21, 5.0),
        "x_ticks": (0.3, 0.5, 1.0, 2.0, 5.0),
        "y_range": (1.0e-7, 1.0),
        "legend_loc": "lower right",
        "run_label_tag": "PBC BC10 (fermiophilic ALP)",
        "figure": "bc10_grendel_talk_zoom",
        "x_clip": (0.3, 4.2),   # GKOZ figure frame: m_a in [0.3, 4.2] GeV
    },
}


def insert_threshold_tips(rows: list[dict[str, object]], threshold: float | None = None) -> list[dict[str, object]]:
    """Close every finite sensitivity transition at ``peak_N == threshold``.

    The two coupling roots merge at the maximum of the yield scan. The
    crossing is interpolated in log yield and log coupling whenever adjacent
    supported rows change sensitivity state. Rows with a non-finite peak (the
    explicitly unsupported eta and eta-prime windows) are never bridged. A
    curve solved at a different signal threshold closes its tips at that
    threshold, not at 3.
    """
    threshold = SIGNAL_THRESHOLD if threshold is None else float(threshold)
    rows = sorted((row.copy() for row in rows), key=lambda row: float(row["mass_GeV"]))
    if len(rows) < 2:
        return rows
    output: list[dict[str, object]] = []
    for left, right in zip(rows[:-1], rows[1:]):
        output.append(left)
        if bool(left["has_sensitivity"]) == bool(right["has_sensitivity"]):
            continue
        if not (np.isfinite(left["peak_N"]) and np.isfinite(right["peak_N"])
                and np.isfinite(left["peak_c"]) and np.isfinite(right["peak_c"])
                and float(left["peak_N"]) > 0 and float(right["peak_N"]) > 0
                and float(left["peak_c"]) > 0 and float(right["peak_c"]) > 0
                and float(right["mass_GeV"]) > float(left["mass_GeV"])):
            continue
        log_left = np.log10(float(left["peak_N"]))
        log_right = np.log10(float(right["peak_N"]))
        log_threshold = np.log10(threshold)
        if log_left == log_right:
            continue
        frac = float((log_threshold - log_left) / (log_right - log_left))
        if not 0.0 < frac < 1.0:
            continue
        tip_c = 10.0 ** (np.log10(float(left["peak_c"]))
                         + frac * (np.log10(float(right["peak_c"])) - np.log10(float(left["peak_c"]))))
        tip = (left if bool(left["has_sensitivity"]) else right).copy()
        tip["mass_GeV"] = float(left["mass_GeV"] + frac * (float(right["mass_GeV"]) - float(left["mass_GeV"])))
        tip["c_min"] = tip_c
        tip["c_max"] = tip_c
        tip["peak_N"] = threshold
        tip["peak_c"] = tip_c
        tip["has_sensitivity"] = True
        output.append(tip)
    output.append(rows[-1])
    return sorted(output, key=lambda row: float(row["mass_GeV"]))


def sensitive_segments(rows: list[dict[str, object]], min_width_GeV: float = 0.0) -> list[list[dict[str, object]]]:
    """Sensitive runs without bridging explicit insensitive rows.

    ``min_width_GeV`` drops islands narrower than that mass width: a band a
    couple of grid steps wide whose peak yield sits a few percent above the
    threshold is a marginal crossing that renders as a detached speck.
    Dropped islands are reported on stderr, never removed silently.
    """
    segments: list[list[dict[str, object]]] = []
    current: list[dict[str, object]] = []
    for row in sorted(rows, key=lambda item: float(item["mass_GeV"])):
        if row["has_sensitivity"]:
            current.append(row)
        elif current:
            segments.append(current)
            current = []
    if current:
        segments.append(current)
    if min_width_GeV <= 0.0:
        return segments
    kept: list[list[dict[str, object]]] = []
    for segment in segments:
        masses = [float(row["mass_GeV"]) for row in segment]
        width = max(masses) - min(masses)
        if width < min_width_GeV:
            peaks = [float(row["peak_N"]) for row in segment if not math.isnan(float(row["peak_N"]))]
            print(f"  dropped island {min(masses):.4f}-{max(masses):.4f} GeV (width {width:.4f} < {min_width_GeV:g} GeV"
                  + (f", peak_N max {max(peaks):.3f}" if peaks else "") + ")", file=sys.stderr)
            continue
        kept.append(segment)
    return kept


def load_island_rows(path: Path, config: dict[str, object], *, with_threshold_tips: bool = True,
                     threshold: float | None = None) -> list[dict[str, object]]:
    """The island rows of a sensitivity CSV for this benchmark.

    Open-top rows (the BC4 e e -only points below 2 m_mu) have an upper edge
    above the scanned coupling ceiling; they are drawn up to the ceiling.
    """
    prefix = str(config["prefix"])
    rows: list[dict[str, object]] = []
    with Path(path).open(newline="") as fh:
        for raw in csv.DictReader(fh):
            c_max = parse_csv_float(raw.get(f"{prefix}_max"))
            if math.isnan(c_max) and parse_csv_bool(raw.get(f"{prefix}_max_open")):
                c_max = float(config.get("scan_ceiling", 1.0e-2))
            rows.append({
                "mass_GeV": parse_csv_float(raw.get("mass_GeV")),
                "c_min": parse_csv_float(raw.get(f"{prefix}_min")),
                "c_max": c_max,
                "peak_N": parse_csv_float(raw.get("peak_N")),
                "peak_c": parse_csv_float(raw.get(f"peak_{prefix}")),
                "has_sensitivity": parse_csv_bool(raw.get("has_sensitivity")),
                "c_min_open": parse_csv_bool(raw.get(f"{prefix}_min_open")),
                "c_max_open": parse_csv_bool(raw.get(f"{prefix}_max_open")),
            })
    m_min = config.get("grendel_m_min")
    if m_min is not None:
        rows = [row for row in rows if row["mass_GeV"] >= float(m_min)]
    return insert_threshold_tips(rows, threshold=threshold) if with_threshold_tips else rows


def load_variation_envelope(path: Path, config: dict[str, object],
                            central_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """Load and validate a non-statistical, one-source-at-a-time envelope
    against the central rows it is drawn around."""
    prefix = str(config["prefix"])
    rows: list[dict[str, object]] = []
    with Path(path).open(newline="") as fh:
        for raw in csv.DictReader(fh):
            definition = str(raw.get("envelope_definition", "")).strip()
            if definition != "single_source_variation_envelope":
                raise ValueError(f"{path}: unsupported envelope_definition {definition!r}")
            row: dict[str, object] = {"mass_GeV": parse_csv_float(raw.get("mass_GeV")),
                                      "has_sensitivity": parse_csv_bool(raw.get("has_sensitivity")),
                                      "break_before": False}
            for edge in ("min", "max"):
                stem = f"{prefix}_{edge}"
                row[f"{edge}_central"] = parse_csv_float(raw.get(f"{stem}_central"))
                row[f"{edge}_lo"] = parse_csv_float(raw.get(f"{stem}_envelope_lo"))
                row[f"{edge}_hi"] = parse_csv_float(raw.get(f"{stem}_envelope_hi"))
                row[f"{edge}_open"] = parse_csv_bool(raw.get(f"{stem}_open"))
                # BC4 and BC10 name the same diagnostics slightly differently. A
                # halo is only meaningful where every physical variation supplies
                # a closed boundary and the same-physics repeats are smaller
                # than the physical envelope.
                row[f"{edge}_variation_open"] = (parse_csv_bool(raw.get(f"{stem}_envelope_open"))
                                                 or parse_csv_bool(raw.get(f"{stem}_any_variation_open")))
                row[f"{edge}_variation_missing"] = parse_csv_bool(raw.get(f"{stem}_variation_missing"))
                row[f"{edge}_numerically_unresolved"] = (
                    parse_csv_bool(raw.get(f"{stem}_numerical_repeat_not_subdominant"))
                    or parse_csv_bool(raw.get(f"{stem}_repeat_not_subdominant")))
            rows.append(row)
    # The envelope is validated mass by mass against the central rows, which
    # ``load_island_rows`` has already cut at ``grendel_m_min``; cut the
    # envelope the same way.
    m_min = config.get("grendel_m_min")
    if m_min is not None:
        rows = [row for row in rows if row["mass_GeV"] >= float(m_min)]
    central_by_mass = {round(float(r["mass_GeV"]), 9): r for r in central_rows if np.isfinite(r["mass_GeV"])}
    if len(rows) > len(central_by_mass):
        raise ValueError(f"{path}: envelope has more masses ({len(rows)}) than the central curve "
                         f"({len(central_by_mass)})")
    central_sorted = sorted(central_by_mass.values(), key=lambda item: float(item["mass_GeV"]))
    previous_mass = None
    for row in sorted(rows, key=lambda item: float(item["mass_GeV"])):
        central = central_by_mass.get(round(float(row["mass_GeV"]), 9))
        if central is None:
            raise ValueError(f"{path}: envelope mass {row['mass_GeV']} is not central")
        if bool(row["has_sensitivity"]) != bool(central["has_sensitivity"]):
            raise ValueError(f"{path}: sensitivity state differs at mass {row['mass_GeV']}")
        for edge in ("min", "max"):
            observed = float(row[f"{edge}_central"])
            expected = float(central[f"c_{edge}"])
            observed_open = bool(row[f"{edge}_open"])
            expected_open = bool(central.get(f"c_{edge}_open", False))
            if observed_open != expected_open:
                raise ValueError(f"{path}: {prefix}_{edge}_open differs from the central curve at mass {row['mass_GeV']}")
            values_match = (not np.isfinite(observed) if expected_open or not bool(central["has_sensitivity"])
                            else np.isfinite(expected) and np.isclose(observed, expected, rtol=1.0e-9, atol=0.0))
            if not values_match:
                raise ValueError(f"{path}: {prefix}_{edge}_central differs from the central curve at mass {row['mass_GeV']}")
        if previous_mass is not None:
            row["break_before"] = any(previous_mass < float(item["mass_GeV"]) < float(row["mass_GeV"])
                                      and not bool(item["has_sensitivity"]) for item in central_sorted)
        previous_mass = float(row["mass_GeV"])
    return rows


def envelope_segments(rows: list[dict[str, object]], edge: str) -> list[list[dict[str, object]]]:
    """Split an envelope boundary at every unsupported or insensitive row."""
    segments: list[list[dict[str, object]]] = []
    current: list[dict[str, object]] = []
    for row in sorted(rows, key=lambda item: float(item["mass_GeV"])):
        if bool(row.get("break_before", False)) and current:
            segments.append(current)
            current = []
        valid = (bool(row["has_sensitivity"]) and not bool(row[f"{edge}_open"])
                 and not bool(row[f"{edge}_variation_open"]) and not bool(row[f"{edge}_variation_missing"])
                 and not bool(row[f"{edge}_numerically_unresolved"])
                 and np.isfinite(row[f"{edge}_lo"]) and np.isfinite(row[f"{edge}_hi"])
                 and float(row[f"{edge}_lo"]) > 0 and float(row[f"{edge}_hi"]) > 0)
        if valid:
            current.append(row)
        elif current:
            segments.append(current)
            current = []
    if current:
        segments.append(current)
    return segments


def plot_variation_envelope(ax, config: dict[str, object], rows: list[dict[str, object]]) -> None:
    """A light halo around each central contour boundary."""
    yscale = float(config.get("y_scale", 1.0))
    labelled = False
    for edge in ("min", "max"):
        for segment in envelope_segments(rows, edge):
            x = np.asarray([row["mass_GeV"] for row in segment], dtype=float)
            lo = yscale * np.asarray([row[f"{edge}_lo"] for row in segment], dtype=float)
            hi = yscale * np.asarray([row[f"{edge}_hi"] for row in segment], dtype=float)
            ax.fill_between(x, lo, hi, facecolor="#f28e8b", edgecolor="none", alpha=0.34, zorder=7,
                            label=("GRENDEL single-source\nmodel variations" if not labelled else None))
            labelled = True


def plot_grendel(ax, config: dict[str, object], rows: list[dict[str, object]]) -> None:
    """The filled GRENDEL island(s)."""
    segments = sensitive_segments(rows, float(config.get("min_island_GeV", 0.0)))
    if not segments:
        return
    yscale = float(config.get("y_scale", 1.0))
    for index, segment in enumerate(segments):
        x = np.asarray([row["mass_GeV"] for row in segment], dtype=float)
        lo = yscale * np.asarray([row["c_min"] for row in segment], dtype=float)
        hi = yscale * np.asarray([row["c_max"] for row in segment], dtype=float)
        ax.fill_between(x, lo, hi, color=GRENDEL_RED, alpha=0.18, lw=0, zorder=8,
                        label="GRENDEL\nthis work" if index == 0 else None)
        ax.plot(x, lo, color=GRENDEL_RED, lw=2.2, zorder=9)
        ax.plot(x, hi, color=GRENDEL_RED, lw=2.2, zorder=9)
        if len(segment) == 1:
            ax.plot([x[0], x[0]], [lo[0], hi[0]], color=GRENDEL_RED, lw=2.2, zorder=9)


def plot_grendel_threshold(ax, config: dict[str, object], rows: list[dict[str, object]], label: str,
                           color: str = GRENDEL_DARK_RED, lw: float = 1.7, ls: object = (0, (5, 2))) -> None:
    """A stricter-threshold island as an unfilled dashed contour, drawn on
    top of the filled baseline band it is nested inside."""
    yscale = float(config.get("y_scale", 1.0))
    for index, segment in enumerate(sensitive_segments(rows, float(config.get("min_island_GeV", 0.0)))):
        x = np.asarray([row["mass_GeV"] for row in segment], dtype=float)
        lo = yscale * np.asarray([row["c_min"] for row in segment], dtype=float)
        hi = yscale * np.asarray([row["c_max"] for row in segment], dtype=float)
        ax.plot(x, lo, color=color, lw=lw, ls=ls, zorder=11, label=label if index == 0 else None)
        ax.plot(x, hi, color=color, lw=lw, ls=ls, zorder=11)
        if len(segment) == 1:
            ax.plot([x[0], x[0]], [lo[0], hi[0]], color=color, lw=lw, ls=ls, zorder=11)


def plot_excluded(ax, config: dict[str, object]) -> None:
    segments = load_polygon_segments(Path(config["excluded"]))
    yscale = float(config.get("y_scale", 1.0))
    segments = [np.column_stack([seg[:, 0], yscale * seg[:, 1]]) for seg in segments]
    label = str(config["excluded_label"])
    if config["excluded_closed"]:
        for idx, seg in enumerate(segments):
            ax.fill(seg[:, 0], seg[:, 1], facecolor="0.62", edgecolor="none", alpha=0.38, zorder=0,
                    label=label if idx == 0 else None)
        return
    # Open boundary chains: the excluded region lies above each chain, up to
    # the source figure's frame top (never extrapolated past the published
    # region). Close each chain along that top edge.
    top = max(seg[:, 1].max() for seg in segments)
    clip = config["x_clip"]
    for idx, seg in enumerate(segments):
        x, y = seg[:, 0], seg[:, 1]
        if clip is not None:
            x = np.clip(x, clip[0], clip[1])
        xs = np.concatenate([x, [x[-1], x[0]]])
        ys = np.concatenate([y, [top, top]])
        ax.fill(xs, ys, facecolor="0.62", edgecolor="none", alpha=0.38, zorder=0,
                label=label if idx == 0 else None)


def plot_extra_excluded(ax, config: dict[str, object]) -> None:
    """Additional existing-exclusion regions (gray fills)."""
    yscale = float(config.get("y_scale", 1.0))
    for entry in config.get("extra_excluded", []):
        segments = load_polygon_segments(Path(entry["path"]))
        segments = [np.column_stack([seg[:, 0], yscale * seg[:, 1]]) for seg in segments]
        label = entry.get("label")
        if label is not None:
            label = str(label)
        if entry["mode"] == "closed":
            for idx, seg in enumerate(segments):
                ax.fill(seg[:, 0], seg[:, 1], facecolor="0.62", edgecolor="none", alpha=0.38, zorder=0,
                        label=label if idx == 0 else None)
        elif entry["mode"] == "band":
            lo, hi = segments[0], segments[1]     # excluded band between the lower and upper chain
            xs = np.concatenate([lo[:, 0], hi[::-1, 0]])
            ys = np.concatenate([lo[:, 1], hi[::-1, 1]])
            ax.fill(xs, ys, facecolor="0.62", edgecolor="none", alpha=0.38, zorder=0, label=label)


def plot_projections(ax, config: dict[str, object]) -> None:
    """Competitor sensitivity projections as line curves."""
    yscale = float(config.get("y_scale", 1.0))
    for entry in config.get("projections", []):
        segments = load_polygon_segments(Path(entry["path"]))
        segments = [np.column_stack([seg[:, 0], yscale * seg[:, 1]]) for seg in segments]
        clip = entry.get("x_clip")
        first = True
        for seg in segments:
            x, y = seg[:, 0], seg[:, 1]
            if clip is not None:
                keep = (x >= clip[0]) & (x <= clip[1])
                idx = np.where(keep)[0]
                if idx.size < 2:
                    continue
                splits = np.where(np.diff(idx) > 1)[0]
                for chunk in np.split(idx, splits + 1):       # split on gaps: never bridge a clip
                    if chunk.size < 2:
                        continue
                    ax.plot(x[chunk], y[chunk], color=entry["color"], ls=entry["ls"], lw=1.5, zorder=6,
                            label=str(entry["label"]) if first else None)
                    first = False
            else:
                ax.plot(x, y, color=entry["color"], ls=entry["ls"], lw=1.5, zorder=6,
                        label=str(entry["label"]) if first else None)
                first = False
