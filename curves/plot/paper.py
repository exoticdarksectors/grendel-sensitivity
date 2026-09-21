"""The manuscript figures.

BC4 and BC10: one panel each, the central curves with a compact, directly
documented comparison set (SHiP and CODEX-b for BC4, SHiP for BC10) and the
N_sig >= 10 island dashed inside the N_sig >= 3 band. HNL: a 2x2 grid with
the three single-flavour panels (SHiP, FASER2 and CODEX-b projections) and
the shared legend in the fourth cell.
"""
from __future__ import annotations

import copy
from pathlib import Path

from .data import configure_matplotlib

configure_matplotlib()
import matplotlib.pyplot as plt  # noqa: E402

from . import bc_panel, hnl_panel  # noqa: E402
from .common import add_run_label, run_label  # noqa: E402
from .results import GrendelResults  # noqa: E402

# Second, deliberately conservative signal threshold drawn alongside the
# N_sig >= 3 baseline, matching the Higgs benchmark panels.
SECONDARY_THRESHOLD = 10.0

PAPER_PROJECTION_FILES = {"bc4": {"ship_bc4.dat", "codexb_300fb.dat"}, "bc10": {"ship_alp2.dat"}}
PAPER_EXTRA_EXCLUSION_FILES = {"bc4": set(), "bc10": {"beamdumps_past_alp2.dat"}}

PAPER_PROJECTION_DRAWERS = (hnl_panel.plot_ship_projection, hnl_panel.plot_faser2_projection,
                            hnl_panel.plot_codexb_projection)
PANEL_ORDER = ("100", "010", "001")
REGION_LABELS = {"experimental": "HNLimits experimental bounds",
                 "cosmology": "HNLimits cosmological bounds"}


# ------------------------------------------------------------ BC4 / BC10 --

def paper_config(benchmark: str) -> dict[str, object]:
    """The frozen, central-only comparison profile for the manuscript."""
    config = copy.deepcopy(bc_panel.BENCHMARKS[benchmark])
    config["projections"] = [e for e in config.get("projections", [])
                             if Path(e["path"]).name in PAPER_PROJECTION_FILES[benchmark]]
    config["extra_excluded"] = [e for e in config.get("extra_excluded", [])
                                if Path(e["path"]).name in PAPER_EXTRA_EXCLUSION_FILES[benchmark]]
    if benchmark == "bc10":
        config["secondary_ylabel"] = None    # the manuscript uses the g_Y convention only
    return config


def render_bc_paper(benchmark: str, results: GrendelResults, output_dir: Path,
                    secondary_threshold: float | None = SECONDARY_THRESHOLD, output_suffix: str = "") -> Path:
    config = paper_config(benchmark)
    fig, ax = plt.subplots(figsize=(7.15, 4.7))
    bc_panel.plot_excluded(ax, config)
    bc_panel.plot_extra_excluded(ax, config)
    bc_panel.plot_projections(ax, config)
    central_rows = bc_panel.insert_threshold_tips(
        bc_panel.load_island_rows(results.sensitivity(benchmark), config, with_threshold_tips=False))
    bc_panel.plot_grendel(ax, config, central_rows)
    drew_secondary = False
    if secondary_threshold:
        nsig10 = results.nsig10(benchmark)
        if nsig10 is not None:
            secondary_rows = bc_panel.load_island_rows(nsig10, config, threshold=secondary_threshold)
            if secondary_rows:
                bc_panel.plot_grendel_threshold(
                    ax, config, secondary_rows,
                    label="GRENDEL $N_{\\mathrm{sig}}> " f"{secondary_threshold:g}$")
                drew_secondary = True
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*config["x_range"])
    ax.set_xticks(config["x_ticks"])
    ax.set_xticklabels([f"{tick:g}" for tick in config["x_ticks"]])
    ax.set_ylim(*config["y_range"])
    ax.set_xlabel(config["xlabel"], fontsize=9)
    ax.set_ylabel(config["ylabel"], fontsize=9)
    ax.tick_params(axis="both", which="both", labelsize=7.5)
    # The thresholds are named in the legend; the corner label carries the
    # beam and luminosity only.
    add_run_label(ax, run_label(str(config["run_label_tag"]), threshold_line=""), fontsize=7.5)
    ax.grid(True, which="both", alpha=0.18)
    if config.get("secondary_ylabel"):
        yscale = float(config.get("y_scale", 1.0))
        secax = ax.secondary_yaxis("right", functions=(lambda y: y / yscale, lambda y: y * yscale))
        secax.set_ylabel(str(config["secondary_ylabel"]), fontsize=9)
        secax.tick_params(axis="y", which="both", labelsize=7.5)
    handles, labels = ax.get_legend_handles_labels()

    # Single flat legend row: drop the "\n arXiv:..." provenance sub-line (the
    # sources are cited in the figure caption). With a second threshold drawn
    # the two GRENDEL entries must stay distinguishable.
    def _short(text: str) -> str:
        head = text.split("\n")[0]
        if not head.startswith("GRENDEL"):
            return head
        if head == "GRENDEL":
            return "GRENDEL $N_{\\mathrm{sig}}> 3$" if drew_secondary else "GRENDEL (this work)"
        return head

    labels = [_short(text) for text in labels]
    # Inside the axes, along the bottom: every curve sits in the upper part of
    # the frame, so the lowest decade is empty in both benchmarks.
    ax.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=len(labels), fontsize=7,
              columnspacing=1.1, handlelength=1.6, handletextpad=0.5, borderpad=0.5, framealpha=0.92,
              edgecolor="0.7").set_zorder(20)
    fig.subplots_adjust(left=0.11, right=0.88 if config.get("secondary_ylabel") else 0.97, top=0.97, bottom=0.12)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{benchmark}_grendel_paper{output_suffix}.pdf"
    fig.savefig(output, bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


# ------------------------------------------------------------------- HNL --

def draw_panel(ax, scenario: str, processed: list[dict[str, object]], results: GrendelResults,
               secondary_threshold: float | None = SECONDARY_THRESHOLD) -> None:
    config = hnl_panel.SCENARIOS[scenario]
    seen: set[str] = set()
    for item in processed:
        category = str(item["category"])
        hnl_panel.plot_limit(ax, item, label=REGION_LABELS[category] if category not in seen else None)
        seen.add(category)
    flavor = config["flavor"]
    # With a second threshold drawn, the baseline band must name its own
    # threshold; "this work" alone would leave N_sig >= 3 identified nowhere.
    hnl_panel.plot_grendel(
        ax, hnl_panel.load_flavor_rows(results.sensitivity("hnl"), flavor),
        label="GRENDEL $N_{\\mathrm{sig}}> 3$" if secondary_threshold else "GRENDEL\nthis work")
    if secondary_threshold:
        nsig10 = results.nsig10("hnl")
        if nsig10 is not None:
            hnl_panel.plot_grendel_threshold(
                ax, hnl_panel.load_flavor_rows(nsig10, flavor, threshold=secondary_threshold),
                label="GRENDEL $N_{\\mathrm{sig}}> " f"{secondary_threshold:g}$")
    for draw_projection in PAPER_PROJECTION_DRAWERS:
        draw_projection(ax, flavor)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*hnl_panel.ZOOM_X_RANGE)
    ax.set_ylim(*hnl_panel.Y_RANGE)
    ax.set_xticks(hnl_panel.ZOOM_X_TICKS)
    ax.set_xticklabels([f"{tick:g}" for tick in hnl_panel.ZOOM_X_TICKS])
    ax.set_xlabel(r"$m_N$ (GeV)")
    ax.set_ylabel(config["latex"])
    ax.grid(True, which="both", alpha=0.18)
    ax.tick_params(axis="both", which="major", labelsize=7)


def processed_by_scenario() -> dict[str, list[dict[str, object]]]:
    return {scenario: hnl_panel.processed_curves(scenario, hnl_panel.SCENARIOS[scenario],
                                                 hnl_panel.load_metadata(hnl_panel.SCENARIOS[scenario]["sheet"]))
            for scenario in PANEL_ORDER}


def render_hnl_paper(results: GrendelResults, output_dir: Path, stem: str = "hnlimits_grendel_paper",
                     secondary_threshold: float | None = SECONDARY_THRESHOLD) -> Path:
    processed = processed_by_scenario()
    # 2x2: three benchmark panels plus the legend in the fourth cell, which
    # keeps the panels the same size as a 2x4 layout at 8 pt legend type and
    # a 5.6 in figure height.
    fig = plt.figure(figsize=(7.2, 5.6))
    grid = fig.add_gridspec(2, 2, left=0.085, right=0.985, top=0.975, bottom=0.095, wspace=0.28, hspace=0.34)
    axes = {"100": fig.add_subplot(grid[0, 0]), "010": fig.add_subplot(grid[0, 1]),
            "001": fig.add_subplot(grid[1, 0])}
    legend_ax = fig.add_subplot(grid[1, 1])
    legend_ax.axis("off")
    for scenario in PANEL_ORDER:
        draw_panel(axes[scenario], scenario, processed[scenario], results, secondary_threshold)
    handles, labels, seen_labels = [], [], set()
    for scenario in PANEL_ORDER:
        for handle, label in zip(*axes[scenario].get_legend_handles_labels()):
            if label and label not in seen_labels:
                handles.append(handle)
                labels.append(label.split("\n")[0])   # the arXiv sub-line is in the caption
                seen_labels.add(label)
    # The run conditions ride as the legend title rather than a suptitle:
    # they are a property of the GRENDEL curves alone (SHiP is an SPS beam
    # dump, FASER2 sits in the forward region).
    legend_ax.legend(handles, labels, loc="center", ncol=1, fontsize=8,
                     title=(r"GRENDEL: HL-LHC, $\mathcal{L}=3\,\mathrm{ab}^{-1}$" "\n"
                            r"$pp$, $\sqrt{s}=14$ TeV"),
                     title_fontsize=8, handlelength=2.4, handletextpad=0.7, labelspacing=0.85,
                     borderpad=0.9, framealpha=0.95)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{stem}.pdf"
    fig.savefig(output)
    fig.savefig(output.with_suffix(".png"), dpi=220)
    plt.close(fig)
    return output
