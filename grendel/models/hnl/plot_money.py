"""The HNL sensitivity contours with every uncertainty band drawn on them.

Built on ``plot_exclusion`` (excluded-region fill, open-edge markers, the
scan ceiling, segment-aware drawing). On top of that it layers the bands
the campaigns produce, each read from the results directory when present:

  * lower edge: FONLL production band (orange, ``hnl_band_fonll.csv``),
    Bc normalisation band (teal, ``bc_nuisance_band.csv``), charged-kaon
    transport band (purple, ``kaon_desc_band.csv``);
  * upper edge: decay-model width band (blue, ``decay_model_band.csv``).

Production and HNL -> SM decay physics only: reconstruction, detector
response, background and statistics are idealised (background-free,
N >= 3). A metadata JSON records what was and was not banded.

    python -m grendel.models.hnl.plot_money --results-dir RESULTS [--out-dir DIR]

expects ``RESULTS/hnl/sensitivity.csv`` and the band CSVs under
``RESULTS/hnl/band/``; ``--central-csv`` and the ``--*-band`` options point
at other files.
"""
from __future__ import annotations

import argparse
import json
import shutil
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .plot_exclusion import _plot_single_panel  # noqa: E402

LAB = {"Ue": r"$|U_e|^2$", "Umu": r"$|U_\mu|^2$", "Utau": r"$|U_\tau|^2$"}
BAND_FILES = {
    "fonll": "hnl_band_fonll.csv",
    "decay": "decay_model_band.csv",
    "bc": "bc_nuisance_band.csv",
    "kaon": "kaon_desc_band.csv",
    "breakdown": "channel_breakdown_u2min.csv",
}


def _provenance(have_fonll, have_kaon):
    fonll = ("FONLL production (orange)" if have_fonll
             else "FONLL production variation not on this figure")
    kaon = ("charged-kaon transport d_esc=[1,3] m (purple)" if have_kaon
            else "charged-kaon transport variation not on this figure")
    return (
        f"Theory/model variations: {fonll}; {kaon}; Bc direct normalisation ±40% "
        "(teal, LHCb arXiv:1910.13404; the induced-tau Bc share is negligible); "
        "HNL total-width/lifetime duality δ(m) (blue, cap 20%, conservative against a ~10% post-QCD residual).  "
        "Not banded: production form factors, absolute visible-BR normalisation, neutral kaons, "
        "detailed material/magnetic transport beyond the d_esc proxy, FONLL αs.  "
        "Reconstruction / detector / background / statistics idealised (background-free, N≥3).  "
        "The island closes on refined grid points where peak_N crosses 3 (3.62-3.70 GeV): "
        "~3.62 GeV (Ue/Umu), ~3.69 GeV (Utau) -- no extrapolation.  "
        "Nuisance-induced topology changes are recorded in the tables but not drawn as ordinary ribbons."
    )


def _dex_densify(bm, b_cen, b_lo, b_hi, cm, c_edge):
    """Densify one band edge onto the fine central grid.

    The band anchors are too coarse to plot as absolute edges -- straight
    log-segments between anchors detach from the wiggly dense central line.
    So the band *width in dex* (relative to the band-run central) is
    interpolated and re-applied to the dense plotted central ``c_edge``, so
    the ribbon hugs every wiggle of the central curve and ends where it does.
    """
    ok = (np.isfinite(b_cen) & (b_cen > 0) & np.isfinite(b_lo) & (b_lo > 0)
          & np.isfinite(b_hi) & (b_hi > 0))
    lo = np.full(len(cm), np.nan)
    hi = np.full(len(cm), np.nan)
    # Interpolate each contiguous finite anchor segment independently: never
    # extrapolate a ribbon through a point where a variation creates or
    # destroys the island.
    valid = np.flatnonzero(ok)
    for run in np.split(valid, np.flatnonzero(np.diff(valid) > 1) + 1):
        if len(run) < 2:
            continue
        support = (cm >= bm[run[0]]) & (cm <= bm[run[-1]])
        dex_lo = np.log10(b_cen[run] / b_lo[run])
        dex_hi = np.log10(b_hi[run] / b_cen[run])
        lo[support] = c_edge[support] * 10.0 ** (-np.interp(cm[support], bm[run], dex_lo))
        hi[support] = c_edge[support] * 10.0 ** (+np.interp(cm[support], bm[run], dex_hi))
    return cm, lo, hi


def _ribbon(ax, m, lo, hi, color, label, zorder):
    good = np.isfinite(lo) & np.isfinite(hi) & (lo > 0) & (hi > 0)
    if good.any():
        ax.fill_between(m, np.where(good, lo, np.nan), np.where(good, hi, np.nan),
                        alpha=0.5, color=color, linewidth=0, zorder=zorder, label=label)


def _metadata(central_csv, l_int_fb, p_cut_mev, have_fonll, have_kaon):
    in_band = ["HNL total-width/lifetime duality (decay-model band)",
               "direct Bc normalisation"]
    if have_fonll:
        in_band.insert(0, "FONLL heavy-flavour production (scale/PDF/m_Q)")
    if have_kaon:
        in_band.append("charged-kaon transport escape length d_esc in [1,3] m")
    lim = ["production form factors", "absolute visible-BR normalisation",
           "neutral-kaon contribution",
           "detailed material/magnetic transport beyond the d_esc proxy",
           "FONLL alpha_s (sub-dominant)"]
    if not have_fonll:
        lim.insert(0, "FONLL scale/PDF/m_Q band -- not on this figure")
    if not have_kaon:
        lim.insert(0, "charged-kaon d_esc transport band -- not on this figure")
    band_sources = {
        "fonll_lower_edge": ("hnl_band_fonll.csv (grendel.band.hnl.campaign + combine)"
                             if have_fonll else "absent"),
        "decay_model_upper_edge": "decay_model_band.csv (grendel.band.hnl.decay_model_band)",
        "bc_lower_edge": "bc_nuisance_band.csv (grendel.band.hnl.bc_nuisance, SIGMA_BC_REL_UNCERT=0.40)",
    }
    if have_kaon:
        band_sources["kaon_lower_edge"] = "kaon_desc_band.csv (d_esc=1/1.5/3 m; nominal 1.5 m)"
    return {
        "result": "GRENDEL HNL theory/model variation diagnostic (single-flavor)",
        "hypothesis": {"mixing": "single-flavor", "flavors": ["Ue", "Umu", "Utau"],
                       "observable": "|U_alpha|^2", "nature": "Majorana",
                       "charge_conjugate_counting": "included (NDecayWidth x2/channel; production both charges)",
                       "production_mixing_equals_decay_mixing": True},
        "beam": {"collision": "pp", "sqrt_s_TeV": 14, "L_int_fb": l_int_fb},
        "limit": {"method": "background-free", "criterion": "N_signal >= 3",
                  "track_momentum_cut_MeV": p_cut_mev},
        "scope": {"in_band": in_band,
                  "interpretation": "theory/model variation diagnostic, not a statistical confidence band",
                  "idealized": ["reconstruction", "detector response", "background", "statistics"]},
        "central_csv": str(central_csv),
        "band_sources": band_sources,
        "limitations_not_banded": lim,
        "known_features": ["Umu/Utau low-mass step where N->l pi closes (m_mu+m_pi=0.245 GeV) -- real physics",
                           "island closes on refined grid points (3.62-3.70 GeV) where peak_N crosses 3: "
                           "~3.62 GeV (Ue/Umu), ~3.69 GeV (Utau) -- no extrapolation"],
        "generated": str(date.today()),
        "reproduce": ["python -m grendel.models.hnl.scan --flavor Ue Umu Utau",
                      "python -m grendel.band.hnl.campaign ...; python -m grendel.band.hnl.combine ...",
                      "python -m grendel.band.hnl.decay_model_band --out band/decay_model_band.csv",
                      "python -m grendel.band.hnl.bc_nuisance --out band/bc_nuisance_band.csv",
                      "python -m grendel.models.hnl.channel_breakdown --out band/channel_breakdown_u2min.csv",
                      "python -m grendel.models.hnl.plot_money --results-dir RESULTS"],
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", type=Path, default=None,
                    help="results tree: <dir>/hnl/sensitivity.csv and <dir>/hnl/band/*.csv")
    ap.add_argument("--central-csv", type=Path, default=None)
    for key, name in BAND_FILES.items():
        ap.add_argument(f"--{key}-band" if key != "breakdown" else "--breakdown", type=Path,
                        default=None, help=f"default: <results>/hnl/band/{name}")
    ap.add_argument("--l-int-fb", type=float, default=3000.0)
    ap.add_argument("--p-cut-mev", type=int, default=100)
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="default: <results>/hnl/money_plot")
    a = ap.parse_args(argv)
    if a.results_dir is None and a.central_csv is None:
        ap.error("give --results-dir or --central-csv")

    base = a.results_dir / "hnl" if a.results_dir else a.central_csv.parent
    central_csv = a.central_csv or base / "sensitivity.csv"
    band_dir = base / "band"
    paths = {key: (getattr(a, f"{key}_band") if key != "breakdown" else a.breakdown)
             or band_dir / name for key, name in BAND_FILES.items()}

    cen = pd.read_csv(central_csv)
    have_fonll = paths["fonll"].exists()
    have_kaon = paths["kaon"].exists()
    fonll = pd.read_csv(paths["fonll"]) if have_fonll else None
    dm = pd.read_csv(paths["decay"]) if paths["decay"].exists() else None
    bc = pd.read_csv(paths["bc"]) if paths["bc"].exists() else None
    kb = pd.read_csv(paths["kaon"]) if have_kaon else None

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.6), sharey=True)
    for ax, fl in zip(axes, ["Ue", "Umu", "Utau"]):
        full = cen[cen.flavor == fl].sort_values("mass_GeV")
        sens = full[full.has_sensitivity == True]                            # noqa: E712
        cm = sens["mass_GeV"].to_numpy(float)
        c_min = sens["u2_min"].to_numpy(float)
        c_max = sens["u2_max"].to_numpy(float)
        _plot_single_panel(ax, cen, fl, is_leftmost=(fl == "Ue"), band_df=None)

        if have_fonll:
            bsub = fonll[fonll.flavor == fl].sort_values("mass_GeV")
            bm = bsub["mass_GeV"].to_numpy(float)
            labelled = False
            for col, c_edge in (("u2_min", c_min), ("u2_max", c_max)):
                m, lo, hi = _dex_densify(
                    bm, bsub[f"{col}_central"].to_numpy(float),
                    bsub[f"{col}_band_lo"].to_numpy(float),
                    bsub[f"{col}_band_hi"].to_numpy(float), cm, c_edge)
                _ribbon(ax, m, lo, hi, "orange",
                        None if labelled else "FONLL production variation", zorder=4)
                labelled = True
        if dm is not None:
            dmf = dm[dm.flavor == fl].sort_values("mass_GeV")
            m, lo, hi = _dex_densify(
                dmf["mass_GeV"].to_numpy(float), dmf["u2_max"].to_numpy(float),
                dmf["u2_max_dm_lo"].to_numpy(float), dmf["u2_max_dm_hi"].to_numpy(float),
                cm, c_max)
            _ribbon(ax, m, lo, hi, "steelblue", "decay-model variation (upper)", zorder=4)
        if bc is not None:
            bcf = bc[bc.flavor == fl].sort_values("mass_GeV")
            # Bc-up strengthens the limit (lower boundary); Bc-down weakens it.
            m, lo, hi = _dex_densify(
                bcf["mass_GeV"].to_numpy(float), bcf["u2_min"].to_numpy(float),
                bcf["u2_min_bc_hi"].to_numpy(float), bcf["u2_min_bc_lo"].to_numpy(float),
                cm, c_min)
            _ribbon(ax, m, lo, hi, "teal", "Bc normalisation variation (lower)", zorder=4)
        if kb is not None:
            kbf = kb[kb.flavor == fl].sort_values("mass_GeV")
            if len(kbf):
                m, lo, hi = _dex_densify(
                    kbf["mass_GeV"].to_numpy(float),
                    kbf["u2_min_desc1p5"].to_numpy(float),   # nominal d_esc = 1.5 m
                    kbf["u2_min_desc3p0"].to_numpy(float),   # 3 m strengthens (lower)
                    kbf["u2_min_desc1p0"].to_numpy(float),   # 1 m weakens (upper)
                    cm, c_min)
                _ribbon(ax, m, lo, hi, "purple", "kaon transport d_esc [1,3] m (lower)", zorder=5)

        ax.set_title(f"HNL {LAB[fl]}", fontsize=13)
        if fl in ("Umu", "Utau"):
            ax.annotate(r"$N\to\ell\pi$ closes" + "\n→ weaker limit", xy=(0.25, 3e-6),
                        xytext=(0.33, 5e-5), fontsize=7, color="dimgray",
                        arrowprops=dict(arrowstyle="->", color="dimgray", lw=0.8))
        ax.legend(fontsize=7, loc="lower right")
    fig.suptitle("GRENDEL HNL theory/model variation diagnostics (single-flavor, Majorana, 14 TeV, "
                 f"{a.l_int_fb:.0f} fb$^{{-1}}$, P>{a.p_cut_mev} MeV)", y=1.0, fontsize=13)
    fig.text(0.5, 0.005, _provenance(have_fonll, have_kaon), ha="center", va="bottom",
             fontsize=6.2, style="italic", wrap=True)
    plt.tight_layout(rect=[0, 0.075, 1, 0.97])

    out = a.out_dir or base / "money_plot"
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "money_plot.png", dpi=130, bbox_inches="tight")
    fig.savefig(out / "money_plot.pdf", bbox_inches="tight")
    (out / "run_metadata.json").write_text(json.dumps(
        _metadata(central_csv, a.l_int_fb, a.p_cut_mev, have_fonll, have_kaon), indent=2))
    shutil.copy2(central_csv, out / "hnl_sensitivity_central.csv")
    for key, name in BAND_FILES.items():
        if paths[key].exists():
            shutil.copy2(paths[key], out / name)
    print(f"money plot + tables -> {out}/  "
          f"(FONLL band: {'yes' if have_fonll else 'absent'}; kaon band: {'yes' if have_kaon else 'absent'})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
