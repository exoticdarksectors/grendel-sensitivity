"""Diagnostic plot of the BC4 exclusion island.

Optionally overlays existing bounds and projections from CSVs in
``data/competitors/`` (one per experiment, columns ``m_S_GeV,sin2theta``);
missing curves are skipped with a note. The comparison figures of the paper
are drawn by ``curves.plot`` from the digitised reference curves instead.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PLOT_S2T_MIN, PLOT_S2T_MAX = 1e-12, 1e-4
PLOT_M_MIN, PLOT_M_MAX = 0.1, 5.0
SIGNAL_THRESHOLD = 3.0

_COMP_DIR = Path(__file__).resolve().parent / "data" / "competitors"

# name -> (label, color, style).  Solid = existing bound, dashed = projection.
_COMPETITORS = {
    "CHARM":    ("CHARM",            "0.35",     "-"),
    "LHCb_BKmumu": ("LHCb B$\\to$K$\\mu\\mu$", "darkgreen", "-"),
    "MATHUSLA": ("MATHUSLA",         "tab:blue",   "--"),
    "CODEXb":   ("CODEX-b",          "tab:purple", "--"),
    "ANUBIS":   ("ANUBIS",           "tab:cyan",   "--"),
    "SHiP":     ("SHiP",             "tab:orange", "--"),
}


def _segments(df):
    """Yield contiguous sensitive runs (no bridging insensitive masses)."""
    sens = df["has_sensitivity"].fillna(False).astype(bool)
    grp = sens.ne(sens.shift(fill_value=False)).cumsum()
    for _, seg in df[sens].groupby(grp[sens]):
        yield seg.sort_values("mass_GeV")


def _with_threshold_tips(df):
    """Insert zero-width contour tips at finite ``peak_N == 3`` crossings."""
    rows = df.sort_values("mass_GeV").to_dict("records")
    if len(rows) < 2:
        return df.sort_values("mass_GeV")
    output = []
    for left, right in zip(rows[:-1], rows[1:]):
        output.append(left)
        if bool(left["has_sensitivity"]) == bool(right["has_sensitivity"]):
            continue
        if not all(
            np.isfinite(float(row[field])) and float(row[field]) > 0
            for row in (left, right)
            for field in ("peak_N", "peak_u2")
        ):
            continue
        log_left = np.log10(float(left["peak_N"]))
        log_right = np.log10(float(right["peak_N"]))
        if log_left == log_right:
            continue
        frac = (
            np.log10(SIGNAL_THRESHOLD) - log_left
        ) / (log_right - log_left)
        if not 0.0 < frac < 1.0:
            continue
        coupling = 10.0 ** (
            np.log10(float(left["peak_u2"]))
            + frac
            * (np.log10(float(right["peak_u2"])) - np.log10(float(left["peak_u2"])))
        )
        tip = dict(left if bool(left["has_sensitivity"]) else right)
        tip.update({
            "mass_GeV": float(left["mass_GeV"])
            + frac * (float(right["mass_GeV"]) - float(left["mass_GeV"])),
            "u2_min": coupling,
            "u2_max": coupling,
            "u2_min_open": False,
            "u2_max_open": False,
            "peak_N": SIGNAL_THRESHOLD,
            "peak_u2": coupling,
            "has_sensitivity": True,
        })
        output.append(tip)
    output.append(rows[-1])
    return pd.DataFrame(output, columns=df.columns).sort_values("mass_GeV")


def _overlay_competitors(ax):
    drawn, missing = [], []
    for name, (label, color, style) in _COMPETITORS.items():
        path = _COMP_DIR / f"{name}.csv"
        if not path.exists() or path.stat().st_size == 0:
            missing.append(name)
            continue
        c = pd.read_csv(path)
        ax.plot(c["m_S_GeV"], c["sin2theta"], style, color=color, lw=1.5,
                label=label, zorder=3)
        drawn.append(name)
    if missing:
        print(f"  (competitor curves not found, skipped: {', '.join(missing)};\n"
              f"   add CSVs to {_COMP_DIR})")
    return drawn


def plot_island(island_csv, output_dir=None, basename="bc4_exclusion"):
    island_csv = Path(island_csv)
    df = _with_threshold_tips(pd.read_csv(island_csv))
    output_dir = Path(output_dir) if output_dir else island_csv.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    _overlay_competitors(ax)

    plotted = False
    for seg in _segments(df):
        mass = seg["mass_GeV"].to_numpy(float)
        u2_min = seg["u2_min"].to_numpy(float)
        u2_max = seg["u2_max"].to_numpy(float)
        lo_open = seg.get("u2_min_open", pd.Series(False, index=seg.index)).fillna(False).to_numpy(bool)
        hi_open = seg.get("u2_max_open", pd.Series(False, index=seg.index)).fillna(False).to_numpy(bool)
        lower = np.where(lo_open, PLOT_S2T_MIN, u2_min)
        upper = np.where(hi_open, PLOT_S2T_MAX, u2_max)
        ax.fill_between(mass, lower, upper, color="red", alpha=0.30,
                        label="GRENDEL (3000 fb$^{-1}$)" if not plotted else None,
                        zorder=5)
        ax.plot(mass, np.where(lo_open, np.nan, u2_min), "r-", lw=1.8, zorder=6)
        ax.plot(mass, np.where(hi_open, np.nan, u2_max), "r-", lw=1.8, zorder=6)
        for end in (0, -1):
            if not lo_open[end] and not hi_open[end]:
                ax.plot([mass[end]] * 2, [u2_min[end], u2_max[end]], "r-", lw=1.8, zorder=6)
        plotted = True

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(PLOT_M_MIN, PLOT_M_MAX)
    ax.set_ylim(PLOT_S2T_MIN, PLOT_S2T_MAX)
    ax.set_xlabel(r"$m_S$ [GeV]", fontsize=14)
    ax.set_ylabel(r"$\sin^2\theta$", fontsize=14)
    ax.set_title("BC4 dark scalar — GRENDEL (PX56, CMS)", fontsize=13)
    ax.grid(True, which="both", alpha=0.2, lw=0.5)
    ax.legend(fontsize=9, loc="lower left", ncol=2)

    out_png = output_dir / f"{basename}.png"
    out_pdf = output_dir / f"{basename}.pdf"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_png}\nSaved: {out_pdf}")
    return out_png
