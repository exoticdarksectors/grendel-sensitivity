"""Wide layouts for projection.

The same loaders and drawing helpers as the paper figures, with the full
diagnostic comparison set (every projection in ``bc_panel.BENCHMARKS``,
ANUBIS on the HNL panels), a larger frame and, on request, the
single-source variation envelope around the BC4/BC10 island.

Do not add or remove physics layers here: a talk panel that must show
something the paper panel does not should change the paper renderer, so
the two cannot diverge.
"""
from __future__ import annotations

from pathlib import Path

from .data import configure_matplotlib

configure_matplotlib()
import matplotlib.pyplot as plt  # noqa: E402

from . import bc_panel, hnl_panel  # noqa: E402
from .common import add_run_label, run_label  # noqa: E402
from .paper import PANEL_ORDER, SECONDARY_THRESHOLD, draw_panel, processed_by_scenario  # noqa: E402
from .results import GrendelResults  # noqa: E402

# PBC benchmark number carried by each single-flavour dominance scenario.
BENCHMARK = {"100": "BC6", "010": "BC7", "001": "BC8"}
RUN_CONDITIONS = (r"GRENDEL: HL-LHC, $\mathcal{L}=3\,\mathrm{ab}^{-1}$" "\n"
                  r"$pp$, $\sqrt{s}=14$ TeV")


def render_bc_talk(benchmark: str, results: GrendelResults, output_dir: Path,
                   with_envelope: bool = False) -> Path:
    config = bc_panel.BENCHMARKS[benchmark]
    fig, ax = plt.subplots(figsize=(12.0, 5.8))
    bc_panel.plot_excluded(ax, config)
    bc_panel.plot_extra_excluded(ax, config)
    bc_panel.plot_projections(ax, config)
    canonical_rows = bc_panel.load_island_rows(results.sensitivity(benchmark), config, with_threshold_tips=False)
    if with_envelope:
        envelope = bc_panel.load_variation_envelope(results.envelope(benchmark), config, canonical_rows)
        bc_panel.plot_variation_envelope(ax, config, envelope)
    bc_panel.plot_grendel(ax, config, bc_panel.insert_threshold_tips(canonical_rows))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*config["x_range"])
    ax.set_xticks(config["x_ticks"])
    ax.set_xticklabels([f"{tick:g}" for tick in config["x_ticks"]])
    ax.set_ylim(*config["y_range"])
    ax.set_xlabel(config["xlabel"])
    ax.set_ylabel(config["ylabel"])
    add_run_label(ax, run_label(str(config["run_label_tag"])))
    ax.grid(True, which="both", alpha=0.18)
    if config.get("secondary_ylabel"):
        yscale = float(config.get("y_scale", 1.0))
        secax = ax.secondary_yaxis("right", functions=(lambda y: y / yscale, lambda y: y * yscale))
        secax.set_ylabel(str(config["secondary_ylabel"]))
    ax.legend(loc=str(config.get("legend_loc", "lower right")), fontsize=6.5, title=bc_panel.LEGEND_TITLE,
              title_fontsize=7, columnspacing=0.9, handlelength=1.6, handletextpad=0.45, labelspacing=0.35,
              borderpad=0.45, framealpha=0.92)
    fig.tight_layout()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{config['figure']}{'_uncertainty' if with_envelope else ''}.pdf"
    fig.savefig(output, bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def render_hnl_talk(results: GrendelResults, output_dir: Path,
                    secondary_threshold: float | None = SECONDARY_THRESHOLD) -> list[Path]:
    """One benchmark per figure, legend beside the panel, larger type."""
    processed = processed_by_scenario()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = []
    for scenario in PANEL_ORDER:
        config = hnl_panel.SCENARIOS[scenario]
        fig = plt.figure(figsize=(9.0, 4.25))
        grid = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.46], left=0.085, right=0.995, top=0.975,
                                bottom=0.145, wspace=0.04)
        ax = fig.add_subplot(grid[0, 0])
        legend_ax = fig.add_subplot(grid[0, 1])
        legend_ax.axis("off")
        draw_panel(ax, scenario, processed[scenario], results, secondary_threshold)
        # draw_panel sizes type for a 3.6 in paper panel; this one is 5.5 in
        # wide and will be projected.
        ax.tick_params(axis="both", which="major", labelsize=10)
        ax.xaxis.label.set_size(13)
        ax.yaxis.label.set_size(13)
        ax.set_title(f"{BENCHMARK[scenario]} — {config['description']}", fontsize=12, pad=8)
        handles, labels, seen = [], [], set()
        for handle, label in zip(*ax.get_legend_handles_labels()):
            if label and label not in seen:
                handles.append(handle)
                labels.append(label.split("\n")[0])   # the arXiv sub-line goes on the slide
                seen.add(label)
        legend_ax.legend(handles, labels, loc="center", ncol=1, fontsize=10, title=RUN_CONDITIONS,
                         title_fontsize=10, handlelength=2.4, handletextpad=0.7, labelspacing=0.9,
                         borderpad=0.9, framealpha=0.95)
        out = output_dir / f"hnlimits_grendel_talk_{config['flavor']}.png"
        fig.savefig(out, dpi=240)
        plt.close(fig)
        outputs.append(out)
    return outputs
