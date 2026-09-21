"""HNL theory/model variation diagnostics against the HNLimits landscape.

Three-panel figures (one per flavour) with the nominal GRENDEL boundary and
the variation ribbons of ``grendel.band.hnl`` -- the FONLL production
band, the decay-model legs and the B_c normalisation nuisance -- drawn as
dex widths around the dense central edges. They are diagnostics, not
confidence bands, and are not used in the manuscript.

Inputs, under ``--grendel-dir``::

    hnl/sensitivity.csv
    hnl/band/hnl_band_fonll.csv
    hnl/band/decay_model_band.csv
    hnl/band/bc_nuisance_band.csv
"""
from __future__ import annotations

import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd

from .data import HNL_DATA, configure_matplotlib
from .results import HNL_BAND_FILES, GrendelResults

configure_matplotlib()
import matplotlib.pyplot as plt  # noqa: E402

VECTOR = HNL_DATA / "vector"
FLAVORS = ["Ue", "Umu", "Utau"]
LATEX = {"Ue": r"$|U_e|^2$", "Umu": r"$|U_\mu|^2$", "Utau": r"$|U_\tau|^2$"}
SCENARIO = {"Ue": "(1,0,0)", "Umu": "(0,1,0)", "Utau": "(0,0,1)"}

PROJ_COLOR = {"SHiP": "#DA70D6", "SHiP_TP_optimistic": "#FF1493", "SHiP_TP_conservative": "#C71585",
              "FASER2": "#006400", "MATHUSLA": "#FF8C00", "DUNE": "#1E90FF", "FCCee": "#4682B4",
              "ANUBIS": "#4B0082", "CODEX-b": "#008B8B", "FLArE": "#FF0000"}
PROJ_STYLE = {"SHiP_TP_optimistic": ":", "SHiP_TP_conservative": ":", "FCCee": ":"}
GREN_COLOR = "#E60000"
OPEN_UPPER = 1.0e-1


def break_jumps(x, y, max_decades=1.5):
    """Insert NaN breaks where consecutive points jump by more than
    ``max_decades`` so a polyline does not bridge distinct branches."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    good = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x = np.where(good, x, np.nan)
    y = np.where(good, y, np.nan)
    if len(x) < 2:
        return x, y
    ly = np.log10(np.maximum(y, 1e-300))
    finite_pair = np.isfinite(ly[:-1]) & np.isfinite(ly[1:])
    big = np.where(finite_pair & (np.abs(np.diff(ly)) > max_decades))[0]
    if not len(big):
        return x, y
    out_x, out_y, cuts = [], [], set(big.tolist())
    for i in range(len(x)):
        out_x.append(x[i])
        out_y.append(y[i])
        if i in cuts:
            out_x.append(np.nan)
            out_y.append(np.nan)
    return np.array(out_x), np.array(out_y)


def load_two_section(path):
    """A ``bottom``/``top`` (or ``lower``/``upper``) two-section curve file."""
    if not os.path.exists(path):
        return None, None
    sec = "low"
    low, up = [], []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("# upper") or s.startswith("# top"):
                sec = "up"
                continue
            if s.startswith("# lower") or s.startswith("# bottom"):
                sec = "low"
                continue
            if s.startswith("#"):
                continue
            if not s:
                sec = "up"
                continue
            m, u2 = s.split()
            (low if sec == "low" else up).append((float(m), float(u2)))
    return (np.array(low) if low else None), (np.array(up) if up else None)


def load_island(path):
    """A three-column ``mass u2_min u2_max`` island file."""
    if not os.path.exists(path):
        return None
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            m, lo, hi = line.split()
            rows.append((float(m), float(lo), float(hi)))
    return np.array(rows)


def load_bc7_excluded():
    path = VECTOR / "currently_excluded.dat"
    if not path.exists():
        return []
    subs, cur = [], []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                if cur:
                    subs.append(cur)
                    cur = []
                continue
            if s.startswith("#"):
                continue
            m, u2 = s.split()
            cur.append((float(m), float(u2)))
    if cur:
        subs.append(cur)
    return subs


def excluded_union(flavor, m_grid, category="all"):
    """Lower envelope of the raw HNLimits exclusions of one flavour."""
    files = sorted(glob.glob(str(HNL_DATA / f"hnlimits_{flavor}" / "*.dat")))
    floor = np.full_like(m_grid, np.inf)
    for f in files:
        is_cosmo = os.path.basename(f).startswith("cosmo_")
        if category == "cosmology" and not is_cosmo:
            continue
        if category == "experimental" and is_cosmo:
            continue
        bot, _ = load_two_section(f)
        if bot is None:
            continue
        bm = bot[:, 0]
        bu = bot[:, 1]
        keep = bu < 9e-3
        if keep.sum() < 2:
            continue
        inside = (m_grid >= bm[keep].min()) & (m_grid <= bm[keep].max())
        if not inside.any():
            continue
        lbm = np.log10(bm[keep])
        lbu = np.log10(bu[keep])
        order = np.argsort(lbm)
        interp = np.interp(np.log10(m_grid[inside]), lbm[order], lbu[order])
        floor[inside] = np.minimum(floor[inside], 10**interp)
    floor[np.isinf(floor)] = 1.0
    return floor


def positive(values):
    arr = np.asarray(values, dtype=float)
    return np.where(np.isfinite(arr) & (arr > 0), arr, np.nan)


def display_upper(values, open_flags, default=OPEN_UPPER, open_takes_priority=False):
    arr = positive(values)
    open_flags = np.asarray(open_flags, dtype=bool)
    if open_takes_priority:
        return np.where(open_flags, default, arr)
    return np.where(np.isfinite(arr), arr, np.where(open_flags, default, np.nan))


def dex_densify(bm, b_cen, b_lo, b_hi, cm, c_edge):
    """Apply the band's dex widths to the dense central edge without
    bridging unsupported anchors."""
    bm = np.asarray(bm, dtype=float)
    b_cen = np.asarray(b_cen, dtype=float)
    b_lo = np.asarray(b_lo, dtype=float)
    b_hi = np.asarray(b_hi, dtype=float)
    cm = np.asarray(cm, dtype=float)
    c_edge = np.asarray(c_edge, dtype=float)
    valid_band = (np.isfinite(bm) & np.isfinite(b_cen) & (b_cen > 0) & np.isfinite(b_lo) & (b_lo > 0)
                  & np.isfinite(b_hi) & (b_hi > 0))
    valid_central = np.isfinite(cm) & np.isfinite(c_edge) & (c_edge > 0)
    lo = np.full(len(cm), np.nan)
    hi = np.full(len(cm), np.nan)
    valid = np.flatnonzero(valid_band)
    for run in np.split(valid, np.flatnonzero(np.diff(valid) > 1) + 1):
        if len(run) < 2:
            continue
        order = run[np.argsort(bm[run])]
        support = valid_central & (cm >= bm[order[0]]) & (cm <= bm[order[-1]])
        dex_lo = np.log10(b_cen[order] / b_lo[order])
        dex_hi = np.log10(b_hi[order] / b_cen[order])
        lo[support] = c_edge[support] * 10.0 ** (-np.interp(cm[support], bm[order], dex_lo))
        hi[support] = c_edge[support] * 10.0 ** (np.interp(cm[support], bm[order], dex_hi))
    return cm, lo, hi


def bc_vertical_edges(frame):
    """Lower/upper plotted edges of the B_c nuisance columns. B_c-up
    strengthens the limit and therefore moves the lower boundary down."""
    return (frame["u2_min_bc_hi"].to_numpy(dtype=float), frame["u2_min_bc_lo"].to_numpy(dtype=float))


def ribbon(ax, m, lo, hi, color, label, zorder=9, alpha=0.35):
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    lower = np.minimum(lo, hi)
    upper = np.maximum(lo, hi)
    good = np.isfinite(m) & np.isfinite(lower) & (lower > 0) & np.isfinite(upper) & (upper > 0)
    if good.any():
        ax.fill_between(m, np.where(good, lower, np.nan), np.where(good, upper, np.nan), color=color,
                        alpha=alpha, lw=0, zorder=zorder, label=label)


def setup_axes(ax, title, flavor, xlim=(0.1, 100), ylim=(1e-12, 1e-1)):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"$m_N$ (GeV)")
    ax.set_ylabel(LATEX[flavor])
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.2)


class Inputs:
    def __init__(self, results: GrendelResults):
        self.central = pd.read_csv(results.sensitivity("hnl"))
        self.fonll = pd.read_csv(results.band("hnl", HNL_BAND_FILES["fonll"]))
        self.decay = pd.read_csv(results.band("hnl", HNL_BAND_FILES["decay"]))
        self.bc = pd.read_csv(results.band("hnl", HNL_BAND_FILES["bc"]))
        self.bc7_excluded = load_bc7_excluded()


def plot_reference_context(ax, inputs, flavor, include_projections=True, include_bc7_excluded=True,
                           split_hnlimits_exclusions=False, label_excluded=False):
    m_grid = np.logspace(-1, 2, 400)
    if split_hnlimits_exclusions:
        for category, color, label, zorder in (("experimental", "0.62", "experimental exclusions", 0),
                                               ("cosmology", "#9EC9E2", "cosmology", 1)):
            ax.fill_between(m_grid, excluded_union(flavor, m_grid, category=category), 1e-2, color=color,
                            alpha=0.48, lw=0, zorder=zorder, label=label)
    else:
        ax.fill_between(m_grid, excluded_union(flavor, m_grid), 1e-2, color="0.6", alpha=0.45, lw=0, zorder=0,
                        label="current exclusions" if label_excluded else None)
    if include_bc7_excluded and flavor == "Umu":
        for sub in inputs.bc7_excluded:
            ax.fill([p[0] for p in sub], [p[1] for p in sub], color="0.55", alpha=0.30, lw=0, zorder=0)
    if not include_projections:
        return
    if flavor == "Umu":
        for name in ["ANUBIS", "CODEX-b", "FLArE", "FASER2", "SHiP"]:
            v = load_island(VECTOR / f"{name}_Umu.dat")
            if v is None:
                continue
            c = PROJ_COLOR[name]
            lm, ly = break_jumps(v[:, 0], v[:, 1])
            um, uy = break_jumps(v[:, 0], v[:, 2])
            ax.plot(lm, ly, "-", color=c, lw=1.2, label=f"{name} (BC7 PDF)")
            ax.plot(um, uy, "-", color=c, lw=1.2)
    src = HNL_DATA / f"hnlimits_projections_{flavor}"
    if src.is_dir():
        for f in sorted(glob.glob(str(src / "*.dat"))):
            name = os.path.basename(f).rsplit("_", 1)[0]
            low, up = load_two_section(f)
            c = PROJ_COLOR.get(name, "#444444")
            style = PROJ_STYLE.get(name, "--")
            if low is not None:
                lm, ly = break_jumps(low[:, 0], low[:, 1])
                ax.plot(lm, ly, style, color=c, lw=0.9, alpha=0.75, label=f"{name} (HNLimits)")
            if up is not None:
                um, uy = break_jumps(up[:, 0], up[:, 1])
                ax.plot(um, uy, style, color=c, lw=0.9, alpha=0.75)


def plot_grendel_band(ax, inputs, flavor):
    central = inputs.central
    g = central[(central["flavor"] == flavor) & central["has_sensitivity"]].sort_values("mass_GeV")
    if g.empty:
        return
    x = positive(g["mass_GeV"])
    lo_c = positive(g["u2_min"])
    max_open = np.asarray(g["u2_max_open"], dtype=bool)
    hi_c_fill = display_upper(g["u2_max"], max_open, open_takes_priority=True)
    hi_c_line = np.where(max_open, np.nan, positive(g["u2_max"]))
    ax.fill_between(x, lo_c, hi_c_fill, color=GREN_COLOR, alpha=0.08, lw=0, zorder=8)

    bsub = inputs.fonll[inputs.fonll["flavor"] == flavor].sort_values("mass_GeV")
    if not bsub.empty:
        bm = bsub["mass_GeV"].to_numpy(dtype=float)
        for i, (col, c_edge) in enumerate((("u2_min", lo_c), ("u2_max", hi_c_line))):
            m, lo, hi = dex_densify(bm, bsub[f"{col}_central"].to_numpy(dtype=float),
                                    bsub[f"{col}_band_lo"].to_numpy(dtype=float),
                                    bsub[f"{col}_band_hi"].to_numpy(dtype=float), x, c_edge)
            ribbon(ax, m, lo, hi, "#F2A000", "FONLL production variation" if i == 0 else None,
                   zorder=9, alpha=0.34)

    dsub = inputs.decay[inputs.decay["flavor"] == flavor].sort_values("mass_GeV")
    if not dsub.empty:
        m, lo, hi = dex_densify(dsub["mass_GeV"].to_numpy(dtype=float), dsub["u2_max"].to_numpy(dtype=float),
                                dsub["u2_max_dm_lo"].to_numpy(dtype=float),
                                dsub["u2_max_dm_hi"].to_numpy(dtype=float), x, hi_c_line)
        ribbon(ax, m, lo, hi, "steelblue", "decay-model variation", zorder=9, alpha=0.34)

    bcsub = inputs.bc[inputs.bc["flavor"] == flavor].sort_values("mass_GeV")
    if not bcsub.empty:
        bc_lower, bc_upper = bc_vertical_edges(bcsub)
        m, lo, hi = dex_densify(bcsub["mass_GeV"].to_numpy(dtype=float), bcsub["u2_min"].to_numpy(dtype=float),
                                bc_lower, bc_upper, x, lo_c)
        ribbon(ax, m, lo, hi, "teal", "Bc normalization variation", zorder=9, alpha=0.34)

    lm, ly = break_jumps(x, lo_c)
    um, uy = break_jumps(x, hi_c_line)
    ax.plot(lm, ly, "-", color=GREN_COLOR, lw=2.2, label="GRENDEL nominal boundary", zorder=10)
    ax.plot(um, uy, "-", color=GREN_COLOR, lw=2.2, zorder=10)


def plot_panel(ax, inputs, flavor, zoom=False, include_projections=True, include_bc7_excluded=True,
               split_hnlimits_exclusions=False, label_excluded=False):
    xlim = (0.25, 4.5) if zoom else (0.1, 100)
    ylim = (1e-10, 1e-1) if zoom else (1e-12, 1e-1)
    setup_axes(ax, f"HNL {LATEX[flavor]} scenario {SCENARIO[flavor]}", flavor, xlim=xlim, ylim=ylim)
    plot_reference_context(ax, inputs, flavor, include_projections=include_projections,
                           include_bc7_excluded=include_bc7_excluded,
                           split_hnlimits_exclusions=split_hnlimits_exclusions, label_excluded=label_excluded)
    plot_grendel_band(ax, inputs, flavor)
    ax.legend(loc="lower left", fontsize=7, ncol=2)


def render_diagnostics(results: GrendelResults, output_dir: Path, benchmarks=("hnl",)) -> list[Path]:
    if "hnl" not in benchmarks:
        return []
    inputs = Inputs(results)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = []

    fig, axes = plt.subplots(1, 3, figsize=(20, 7))
    for ax, flavor in zip(axes, FLAVORS):
        plot_panel(ax, inputs, flavor, zoom=False)
    fig.suptitle("GRENDEL nominal reach with theory/model variation diagnostics vs HNLimits "
                 "(BC7 PDF only for Umu)", y=1.01)
    fig.tight_layout()
    outputs.append(output_dir / "grendel_band_three_flavor.pdf")
    fig.savefig(outputs[-1], bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(20, 7))
    for ax, flavor in zip(axes, FLAVORS):
        plot_panel(ax, inputs, flavor, zoom=True)
        ax.set_title(ax.get_title() + " - zoom")
    fig.tight_layout()
    outputs.append(output_dir / "grendel_band_three_flavor_zoom.pdf")
    fig.savefig(outputs[-1], bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(20, 7))
    for ax, flavor in zip(axes, FLAVORS):
        plot_panel(ax, inputs, flavor, zoom=True, include_projections=False, include_bc7_excluded=False,
                   split_hnlimits_exclusions=True, label_excluded=True)
        ax.set_title(ax.get_title() + " - zoom")
    fig.tight_layout()
    outputs.append(output_dir / "grendel_band_three_flavor_zoom_clean.pdf")
    fig.savefig(outputs[-1], bbox_inches="tight")
    plt.close(fig)

    print("\n--- GRENDEL nominal lower edge by flavor ---")
    print(f'{"m [GeV]":>8s}  {"Ue central":>12s}  {"Umu central":>12s}  {"Utau central":>12s}')
    for m in (0.305, 0.5, 1.0, 1.5, 2.0, 3.0):
        cells = []
        for flavor in FLAVORS:
            g = inputs.central[(inputs.central["flavor"] == flavor) & inputs.central["has_sensitivity"]]
            if len(g) and abs(g["mass_GeV"] - m).min() < max(m * 0.12, 0.02):
                row = g.loc[(g["mass_GeV"] - m).abs().idxmin()]
                cells.append(f'{row["u2_min"]:.2e}')
            else:
                cells.append("--")
        print(f"  {m:>6.3f}  " + "  ".join(f"{c:>12s}" for c in cells))
    print(f"\ncentral mass points: {inputs.central['mass_GeV'].nunique()}")
    return outputs
