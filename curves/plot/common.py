"""Helpers shared by the HNL and BC4/BC10 panels."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

SIGNAL_THRESHOLD = 3.0
GRENDEL_RED = "#d62728"
GRENDEL_DARK_RED = "#7f0d14"


def parse_csv_float(value: object) -> float:
    if value is None:
        return math.nan
    text = str(value).strip()
    return float(text) if text else math.nan


def parse_csv_bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def load_polygon_segments(path: Path) -> list[np.ndarray]:
    """Blank-line-separated ``x y`` segments of a digitised curve file."""
    segments: list[np.ndarray] = []
    rows: list[tuple[float, float]] = []
    for raw_line in Path(path).read_text().splitlines():
        line = raw_line.strip()
        if not line:
            if rows:
                segments.append(np.asarray(rows, dtype=float))
                rows = []
            continue
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            rows.append((float(parts[0]), float(parts[1])))
    if rows:
        segments.append(np.asarray(rows, dtype=float))
    return segments


def run_label(tag: str, threshold_line: str | None = None) -> str:
    """Run-condition box. ``threshold_line`` overrides the last line; an
    empty string omits it, for figures that name their thresholds in the
    legend instead."""
    if threshold_line is None:
        threshold_line = r"$N_{\mathrm{sig}}\geq 3$ (bkg-free)"
    # The run conditions belong to GRENDEL alone. Naming it prevents the box
    # from reading as a property of every curve in the frame: CODEX-b assumes
    # 300/fb, and SHiP is an SPS beam dump rather than pp at 14 TeV.
    lines = [tag,
             r"GRENDEL: HL-LHC, $\mathcal{L}=3\ \mathrm{ab}^{-1}$, "
             r"$pp$ $\sqrt{s}=14\ \mathrm{TeV}$"]
    if threshold_line:
        lines.append(threshold_line)
    return "\n".join(lines)


def add_run_label(ax, text: str, fontsize: float = 11) -> None:
    ax.text(0.98, 0.98, text, transform=ax.transAxes, ha="right", va="top",
            fontsize=fontsize, color="black",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72,
                  "boxstyle": "round,pad=0.25"},
            zorder=20)
