"""Acceptance-weighted production-channel breakdown at the exclusion edges.

Which production channel dominates the *accepted* signal at each boundary,
per flavor and mass. This is not the peak-weighted view, which misleads:
the tau-mixing peak yield is W/Z-dominated, but the lower edge -- the one
production-normalisation nuisances actually move -- is B + Bc.

For each (flavor, mass) the combined scan's ``u2_min / peak_u2 / u2_max``
are taken and, channel by channel, that channel's four-vectors are run
through the same geometry and decay Monte Carlo, evaluating
``N_ch(U^2) = sum_ev weight * P_decay(U^2)`` at those couplings. The
fraction ``N_ch / sum_ch N_ch`` is the acceptance-weighted channel share at
that boundary. Hits are capped (``--max-hit-events``): fractions are
relative, so the cap is unbiased and cheap.

Where the HNL visible fraction is tiny -- muon mixing below the
``N -> mu pi`` threshold (m_N < 0.245 GeV), where only the soft
``N -> mu e nu`` survives -- almost no decay samples pass, per-channel yields
are Monte Carlo noise and the fractions are unreliable. The map is robust
above ~0.3 GeV.

    python -m grendel.models.hnl.channel_breakdown --central-csv DIR/sensitivity.csv --out breakdown.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from ...constants import L_INT_PB
from ...geometry.raycast import compute_geometry, directions_from_eta_phi, get_mesh
from ...io.paths import ModelPaths
from ...io.vectors import format_mass_for_filename, load_combined_csv
from ...reco.acceptance import build_event_mc, scan_u2
from .spec import load_template_bundle, seed_for, select_hit_sample

# low-mass -> high-mass production ordering (for stacked plots)
CHANNELS = ["Kmeson", "Dmeson", "tau", "induced_tau", "Bmeson", "Bbaryon", "Bc", "WZ"]
BOUNDS = ["u2_min", "peak_u2", "u2_max"]
# The default scan excludes 0.2 GeV: below ~0.3 GeV the accepted-hit
# statistics are too sparse for a stable breakdown. It can still be requested
# explicitly with ``--mass 0.2``.
DEFAULT_MASSES = [0.305, 0.5, 0.8, 1.0, 1.4, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6]


def _channel_yields(paths, flavor, mass, mass_label, mesh, templates, u2_targets,
                    max_hit_events, decay_samples):
    """N_ch at each target U^2 for every channel (0 for empty/missing)."""
    ctau = float(templates["ctau_m_u2eq1"])
    u2_grid = np.array(sorted(u2_targets))
    out = {}
    for ch in CHANNELS:
        csv = paths.vectors / flavor / ch / f"mN_{mass_label}.csv"
        out[ch] = dict(zip(u2_grid, np.zeros(len(u2_grid))))
        if not csv.exists() or csv.stat().st_size == 0:
            continue
        data = load_combined_csv(csv, mass)
        if len(data["weight"]) == 0:
            continue
        hits, entry_d, exit_d = compute_geometry(data["eta"], data["phi"], mesh)
        idx_all = np.where(hits & np.isfinite(entry_d) & np.isfinite(exit_d))[0]
        if len(idx_all) == 0:
            continue
        rng = np.random.default_rng(seed_for(flavor, f"{mass_label}_{ch}"))
        idx, w, _ = select_hit_sample(idx_all, data["weight"][idx_all], max_hit_events, rng)
        direction = directions_from_eta_phi(data["eta"][idx], data["phi"][idx])
        p_mag = data["beta_gamma"][idx] * mass
        p4 = np.column_stack([data["gamma"][idx] * mass, p_mag[:, None] * direction])
        d, passed, _ = build_event_mc(p4, direction, entry_d[idx], exit_d[idx],
                                      templates, decay_samples, rng)
        _, N = scan_u2(d, passed, exit_d[idx] - entry_d[idx], w,
                       data["beta_gamma"][idx], ctau, L_INT_PB, u2_grid)
        out[ch] = dict(zip(u2_grid, N))
    return out


def run(paths: ModelPaths, flavors, masses, central_csv, max_hit_events, decay_samples):
    mesh = get_mesh()
    central = pd.read_csv(central_csv)
    central = central[central["has_sensitivity"] == True]                 # noqa: E712
    rows = []
    for flavor in flavors:
        for mass in masses:
            crow = central[(central.flavor == flavor) & (central.mass_GeV == mass)]
            if crow.empty:
                print(f"  {flavor:5s} m={mass}: no central sensitivity -- skip", flush=True)
                continue
            crow = crow.iloc[0]
            targets = {b: float(crow[b]) for b in BOUNDS
                       if np.isfinite(crow[b]) and crow[b] > 0}
            ml = format_mass_for_filename(mass)
            templates = load_template_bundle(paths.templates / flavor / f"templates_{ml}.npz")
            if templates is None:
                print(f"  {flavor:5s} m={mass}: no templates -- skip", flush=True)
                continue
            yields = _channel_yields(paths, flavor, mass, ml, mesh, templates,
                                     targets.values(), max_hit_events, decay_samples)
            for b, u2 in targets.items():
                tot = sum(yields[ch][u2] for ch in CHANNELS)
                for ch in CHANNELS:
                    rows.append({"flavor": flavor, "mass_GeV": mass, "boundary": b,
                                 "u2": u2, "channel": ch,
                                 "N": yields[ch][u2],
                                 "frac": (yields[ch][u2] / tot) if tot > 0 else 0.0})
            if "u2_min" in targets:
                u2m = targets["u2_min"]
                tot = sum(yields[ch][u2m] for ch in CHANNELS)
                if tot > 0:
                    top = sorted(((yields[ch][u2m] / tot, ch) for ch in CHANNELS),
                                 reverse=True)[:2]
                    print(f"  {flavor:5s} m={mass:<4} u2_min top: "
                          + ", ".join(f"{c} {f*100:.0f}%" for f, c in top), flush=True)
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flavor", nargs="+", default=["Ue", "Umu", "Utau"])
    ap.add_argument("--mass", nargs="+", type=float, default=None)
    ap.add_argument("--central-csv", default=None,
                    help="combined-scan sensitivity CSV (default: <analysis>/sensitivity.csv)")
    ap.add_argument("--max-hit-events", type=int, default=2000)
    ap.add_argument("--decay-samples", type=int, default=30)
    ap.add_argument("--vectors-dir", type=Path, default=None)
    ap.add_argument("--templates-dir", type=Path, default=None)
    ap.add_argument("--analysis-dir", type=Path, default=None)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    paths = ModelPaths.resolve("hnl", vectors=args.vectors_dir, templates=args.templates_dir,
                               analysis=args.analysis_dir)
    central_csv = args.central_csv or paths.analysis / "sensitivity.csv"
    df = run(paths, args.flavor, args.mass or DEFAULT_MASSES, central_csv,
             args.max_hit_events, args.decay_samples)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"wrote {args.out} ({len(df)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
