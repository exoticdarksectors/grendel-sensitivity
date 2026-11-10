#!/usr/bin/env python3
"""
hnl/dev/debug_plots.py

Internal production diagnostics from tmp/runs/<tag>/llp_4vectors CSVs.
Dev-only — see hnl/dev/README.md.

  python dev/debug_plots.py --channel-fraction Ue
  python dev/debug_plots.py --channel-fraction all
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HNL_ROOT = Path(__file__).resolve().parent.parent  # hnl/dev/ -> hnl/
sys.path.insert(0, str(HNL_ROOT))

from config_mass_grid import MASS_GRID, format_mass_for_filename
from production.paths import ANALYSIS_DIR, LLP_VECTORS_DIR

OUTPUT_BASE = LLP_VECTORS_DIR
PLOT_DIR = ANALYSIS_DIR / "debug_plots"

FLAVORS = ["Ue", "Umu", "Utau"]
# Same channel folders as combine_channels.py / run_all.py.
FRACTION_CHANNELS = ["Bmeson", "Dmeson", "Bc", "Bbaryon", "tau", "induced_tau", "Kmeson", "WZ"]
CHAN_COLORS = {
    "Bmeson": "C0",
    "Dmeson": "C1",
    "Bc": "C2",
    "Bbaryon": "C7",
    "tau": "C3",
    "induced_tau": "C6",
    "Kmeson": "C4",
    "WZ": "C5",
}
CHAN_LABELS = {
    "Bmeson": "Bmeson",
    "Dmeson": "Dmeson",
    "Bc": "Bc",
    "Bbaryon": r"$\Lambda_b$",
    "tau": r"prompt $\tau$",
    "induced_tau": r"induced $\tau$",
    "Kmeson": r"$K^\pm$",
    "WZ": r"$W/Z$",
}


def _sum_weights(csv_path: Path) -> float:
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return 0.0
    try:
        data = np.loadtxt(csv_path, delimiter=",")
    except ValueError:
        return 0.0
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return float(data[:, 0].sum()) if data.size else 0.0


def channel_weight_table(flavor: str, channels: list[str] | None = None):
    """
    Per-channel sum(weights) on the full MASS_GRID.

    Returns
    -------
    masses : ndarray, shape (N,)
    weights : ndarray, shape (n_channels, N)  — 0 where a point has no CSV / empty file
    """
    channels = channels or FRACTION_CHANNELS
    masses = np.asarray(MASS_GRID, dtype=float)
    weights = np.zeros((len(channels), len(masses)), dtype=float)
    for j, ch in enumerate(channels):
        ch_dir = OUTPUT_BASE / flavor / ch
        for i, m in enumerate(masses):
            csv_path = ch_dir / f"mN_{format_mass_for_filename(m)}.csv"
            weights[j, i] = _sum_weights(csv_path)
    return masses, weights


def active_channels(flavor: str, channels: list[str] | None = None) -> list[str]:
    """Channels with nonzero yield on at least one mass point (skip empty WZ/K, etc.)."""
    channels = channels or FRACTION_CHANNELS
    _, weights = channel_weight_table(flavor, channels)
    return [ch for ch, row in zip(channels, weights) if row.sum() > 0]


def plot_channel_fraction(
    flavor: str,
    out_path: Path | None = None,
    channels: list[str] | None = None,
    *,
    only_active: bool = True,
) -> Path:
    """
    Stacked fraction of total production yield vs m_N on a common mass grid.

    All channels share one sorted x axis; zeros are real (channel closed), not dropped,
    so stackplot never connects disjoint mass points with spurious diagonals.
    """
    channels = list(channels or FRACTION_CHANNELS)
    if only_active:
        channels = active_channels(flavor, channels)
    if not channels:
        raise RuntimeError(f"No active production channels for flavor {flavor}")
    masses, weights = channel_weight_table(flavor, channels)
    total = weights.sum(axis=0)
    positive = total > 0
    if not np.any(positive):
        raise RuntimeError(f"No production weights found for flavor {flavor}")

    m = masses[positive]
    w = weights[:, positive]
    total_pos = total[positive]
    fractions = w / total_pos

    out_path = out_path or (PLOT_DIR / f"channel_fraction_{flavor}.png")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 5))
    colors = [CHAN_COLORS.get(ch, f"C{k}") for k, ch in enumerate(channels)]
    labels = [CHAN_LABELS.get(ch, ch) for ch in channels]
    ax.stackplot(
        m,
        *fractions,
        labels=labels,
        colors=colors,
        alpha=0.85,
        linewidth=0,
    )
    ax.set_xscale("log")
    ax.set_xlim(m.min(), m.max())
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel(r"$m_N$ [GeV]")
    ax.set_ylabel(f"channel fraction of total production ({flavor})")
    ax.set_title(f"Production channel composition vs. mass ({flavor})")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main():
    ap = argparse.ArgumentParser(description="HNL production debug plots")
    ap.add_argument(
        "--channel-fraction",
        nargs="?",
        const="all",
        default=None,
        metavar="FLAVOR",
        help="Plot stacked channel fractions (Ue, Umu, Utau, or all).",
    )
    args = ap.parse_args()

    if args.channel_fraction is None:
        ap.print_help()
        sys.exit(0)

    flavors = FLAVORS if args.channel_fraction == "all" else [args.channel_fraction]
    for flavor in flavors:
        if flavor not in FLAVORS:
            ap.error(f"unknown flavor {flavor!r}; choose from {FLAVORS}")
        path = plot_channel_fraction(flavor)
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
