#!/usr/bin/env python3
"""Draw the BC5 island split by production mode, and at two values of BR(h -> SS).

Reads the BC5 ``sensitivity.csv`` (with the ``mixing_*``, ``hSS_*``, ``BSS_*`` and ``brhss<BR>_*``
columns the scan writes) and draws two panels on the (m_S, sin^2 theta) plane:

  (a) the full island at the nominal BR(h -> SS) and the island of each production mode alone;
  (b) the full island and the h -> SS island at the nominal BR and at the overlay BR.

Only GRENDEL is drawn; other experiments' curves belong to the paper's renderer.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from grendel.constants import L_INT_FB  # noqa: E402
from grendel.io.thresholds import split_threshold  # noqa: E402
from grendel.io.variants import variants_in  # noqa: E402
from grendel.models.scalar.model import M_BPLUS, M_KPLUS  # noqa: E402
from grendel.models.scalar.spec_bc5 import br_tag  # noqa: E402

INK = "#0b0b0b"
INK_SOFT = "#52514e"
GRID = "#d9d8d3"
COLORS = {"hSS": "#2a78d6", "mixing": "#eb6834", "BSS": "#1baf7a"}
LABELS = {"mixing": r"$b \to X_s S$ (mixing only)",
          "hSS": r"$h \to SS$ only",
          "BSS": r"$b \to X_s SS,\ B_s \to SS$ only"}


def island_paths(frame: pd.DataFrame, prefix: str, y_lo: float, y_hi: float):
    """Closed (x, y) outlines of the ``prefix`` island, one per contiguous run of sensitive masses.
    Open edges (the island reaching the scan grid) are closed at ``y_lo`` / ``y_hi``; pass values
    far outside the axes so that the closing segment is clipped away and the island reads as open."""
    head = (prefix + "_") if prefix else ""
    paths = []
    run = []

    def flush():
        if len(run) >= 2:
            m = np.array([r[0] for r in run])
            lo = np.array([r[1] for r in run])
            hi = np.array([r[2] for r in run])
            paths.append((np.concatenate([m, m[::-1], m[:1]]),
                          np.concatenate([lo, hi[::-1], lo[:1]])))
        run.clear()

    for _, row in frame.sort_values("mass_GeV").iterrows():
        if not bool(row.get(f"{head}has_sensitivity", False)):
            flush()
            continue
        lo = row[f"{head}u2_min"]
        hi = row[f"{head}u2_max"]
        if not np.isfinite(lo):
            lo = y_lo if bool(row.get(f"{head}u2_min_open", False)) else np.nan
        if not np.isfinite(hi):
            hi = y_hi if bool(row.get(f"{head}u2_max_open", False)) else np.nan
        if not (np.isfinite(lo) and np.isfinite(hi)):
            flush()
            continue
        run.append((row["mass_GeV"], lo, hi))
    flush()
    return paths


def draw(ax, frame, prefix, *, color, label, fill=None, ls="-", lw=1.8, y_lo=1e-40, y_hi=1e10, z=3):
    first = True
    for x, y in island_paths(frame, prefix, y_lo, y_hi):
        if fill is not None:
            ax.fill(x, y, color=fill, lw=0, zorder=z - 1)
        ax.plot(x, y, color=color, ls=ls, lw=lw, label=label if first else None, zorder=z,
                solid_capstyle="round")
        first = False
    if first:  # nothing drawn: keep the legend entry so the absence is visible
        ax.plot([], [], color=color, ls=ls, lw=lw, label=label + " (none)")


def style(ax, x_lim, y_lim):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*x_lim)
    ax.set_ylim(*y_lim)
    ax.set_xlabel(r"$m_S$ [GeV]", color=INK)
    ax.set_ylabel(r"$\sin^2\theta$", color=INK)
    ax.grid(True, which="major", color=GRID, lw=0.6)
    ax.tick_params(colors=INK_SOFT, which="both", length=3)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_SOFT)
    m_thr = M_BPLUS - M_KPLUS
    ax.axvline(m_thr, color=INK_SOFT, lw=0.8, ls=(0, (1, 3)), zorder=1)
    ax.text(m_thr * 1.08, y_lim[0] * 10 ** (0.025 * np.log10(y_lim[1] / y_lim[0])),
            r"$b \to X_s S$ closed", color=INK_SOFT, fontsize=10.5, va="bottom", ha="left")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sensitivity", type=Path, help="the BC5 sensitivity.csv")
    ap.add_argument("--threshold", type=float, default=None,
                    help="draw the N >= <T> islands (a secondary threshold of the scan; default: primary)")
    ap.add_argument("--overlay", type=float, default=None,
                    help="the overlay BR(h -> SS) for panel (b) (default: the first brhss<BR> variant)")
    ap.add_argument("--out", type=Path, nargs="+", default=None,
                    help="output files (default: bc5_modes.pdf and .png beside the input)")
    ap.add_argument("--xlim", type=float, nargs=2, default=None,
                    help="mass range (default: the scanned masses with a margin)")
    ap.add_argument("--ylim", type=float, nargs=2, default=None,
                    help="sin^2 theta range (default: every drawn island with a decade of margin)")
    args = ap.parse_args(argv)

    frame = pd.read_csv(args.sensitivity)
    if args.threshold is not None:
        frame = split_threshold(frame, args.threshold)
        threshold = args.threshold
    else:
        threshold = None
    br_nominal = float(frame["br_hss"].iloc[0]) if "br_hss" in frame.columns else 0.01
    overlays = sorted({float(v[len("brhss"):]) for v in variants_in(frame)
                       if v.startswith("brhss") and "_" not in v[len("brhss"):]})
    br_overlay = args.overlay if args.overlay is not None else (overlays[0] if overlays else None)

    edge_cols = [c for c in frame.columns if c.endswith(("u2_min", "u2_max"))]
    edges = frame[edge_cols].to_numpy(float)
    edges = edges[np.isfinite(edges) & (edges > 0)]
    x_lim = tuple(args.xlim) if args.xlim else (frame["mass_GeV"].min() / 1.25, frame["mass_GeV"].max() * 1.25)
    if args.ylim:
        y_lim = tuple(args.ylim)
    else:
        lo = 10 ** np.floor(np.log10(edges.min()) - 0.5)
        hi = 10 ** np.ceil(np.log10(edges.max()) + 0.8)   # room for the legend above the islands
        y_lim = (lo, hi)

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # sized for a two-column figure* in a revtex paper: the 11 in wide figure prints at ~0.64 scale
    plt.rcParams.update({"font.size": 13, "axes.labelsize": 15, "legend.fontsize": 11.5,
                         "xtick.labelsize": 12, "ytick.labelsize": 12,
                         "legend.frameon": False, "mathtext.default": "regular"})
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11.0, 4.6), sharey=True)
    kw = {}  # open edges are closed off-canvas (draw's defaults)

    pct = lambda br: f"{100 * br:g}%"  # noqa: E731
    draw(ax_a, frame, "", color=INK, fill="#ececea", lw=2.0,
         label=f"all modes, BR($h \\to SS$) = {pct(br_nominal)}", z=2, **kw)
    for mode in ("hSS", "mixing", "BSS"):
        draw(ax_a, frame, mode, color=COLORS[mode], label=LABELS[mode], **kw)
    ax_a.set_title(f"(a) by production mode, BR($h \\to SS$) = {pct(br_nominal)}", loc="left",
                   fontsize=13, color=INK)

    draw(ax_b, frame, "", color=INK, fill="#ececea", lw=2.0,
         label=f"all modes, BR = {pct(br_nominal)}", z=2, **kw)
    draw(ax_b, frame, "hSS", color=COLORS["hSS"], label=f"$h \\to SS$ only, BR = {pct(br_nominal)}", **kw)
    if br_overlay is not None:
        tag = br_tag(br_overlay)
        draw(ax_b, frame, tag, color=INK, ls="--", lw=1.6, label=f"all modes, BR = {pct(br_overlay)}", **kw)
        draw(ax_b, frame, f"{tag}_hSS", color=COLORS["hSS"], ls="--", lw=1.6,
             label=f"$h \\to SS$ only, BR = {pct(br_overlay)}", **kw)
    draw(ax_b, frame, "mixing", color=COLORS["mixing"], lw=1.2, label=LABELS["mixing"] + ", any BR", **kw)
    ax_b.set_title("(b) scaling with BR($h \\to SS$)", loc="left", fontsize=13, color=INK)

    for ax in (ax_a, ax_b):
        style(ax, x_lim, y_lim)
        ax.legend(loc="upper right", handlelength=2.0, labelspacing=0.35)
    ax_b.set_ylabel("")
    n_text = f"N $\\geq$ {threshold:g}" if threshold is not None else "N $\\geq$ 3"
    fig.text(0.01, 0.01, f"GRENDEL, HL-LHC {L_INT_FB / 1000:g} ab$^{{-1}}$, {n_text} signal events; "
             f"dark scalar, PBC BC5", color=INK_SOFT, fontsize=10.5, ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.03, 1, 1))

    outs = args.out or [args.sensitivity.with_name("bc5_modes.pdf"), args.sensitivity.with_name("bc5_modes.png")]
    for out in outs:
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=200 if out.suffix.lower() == ".png" else None)
        print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
