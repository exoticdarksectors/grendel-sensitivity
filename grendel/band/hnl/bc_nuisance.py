"""Bc production-normalisation nuisance on the HNL lower edge.

The Bc channel owns a large share of the accepted lower-edge yield near the
dome (~28-44% over 1.8-3.6 GeV), and its normalisation is the frozen
``SIGMA_BC_PB``, anchored to the LHCb f(Bc)/f(B) ratio but carrying a real
uncertainty (``SIGMA_BC_REL_UNCERT``, default 0.40). Scaling the Bc channel
by ``(1 +/- delta_Bc)`` moves the lower edge ``u2_min`` wherever Bc
contributes.

Because the signal is additive in channels, this is a linear reweight of
one channel -- no full re-run. For each (flavor, mass) the *combined* and
the *Bc-only* events are scanned over the U^2 grid and

    u2_min_central = boundary(N_tot),
    u2_min_hi      = boundary(N_tot + delta*N_Bc)   (Bc up   -> stronger limit),
    u2_min_lo      = boundary(N_tot - delta*N_Bc)   (Bc down -> weaker  limit).

Scope: this bands the direct Bc channel (Bc -> l N) only. ``SIGMA_BC_PB``
also normalises the Bc -> tau nu source inside ``induced_tau``, but that
source is 0.06% of the tau pool, so co-varying it shifts ``induced_tau`` by
< 0.03% -- negligible against the ~0.2 dex direct-Bc band.

    python -m grendel.band.hnl.bc_nuisance --out bc_nuisance_band.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from ...constants import L_INT_PB, LOG_U2_MAX, LOG_U2_MIN, N_THRESHOLD, N_U2_POINTS
from ...geometry.raycast import compute_geometry, directions_from_eta_phi, get_mesh
from ...io.paths import ModelPaths
from ...io.vectors import format_mass_for_filename, load_combined_csv
from ...models.hnl.spec import load_template_bundle, seed_for, select_hit_sample
from ...production.constants import SIGMA_BC_REL_UNCERT
from ...reco.acceptance import build_event_mc, scan_u2
from ...reco.exclusion import find_exclusion_band

DEFAULT_MASSES = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6]


def _scan_csv(csv, flavor, mass, tag, mesh, templates, u2_grid,
              max_hit_events, decay_samples, seed_salt):
    """N(U^2) for one channel/combined CSV (zeros if empty)."""
    if not csv.exists() or csv.stat().st_size == 0:
        return np.zeros(len(u2_grid))
    data = load_combined_csv(csv, mass)
    if len(data["weight"]) == 0:
        return np.zeros(len(u2_grid))
    hits, entry_d, exit_d = compute_geometry(data["eta"], data["phi"], mesh)
    idx_all = np.where(hits & np.isfinite(entry_d) & np.isfinite(exit_d))[0]
    if len(idx_all) == 0:
        return np.zeros(len(u2_grid))
    rng = np.random.default_rng(seed_for(
        flavor, f"{format_mass_for_filename(mass)}_{tag}", seed_salt))
    idx, w, _ = select_hit_sample(idx_all, data["weight"][idx_all], max_hit_events, rng)
    direction = directions_from_eta_phi(data["eta"][idx], data["phi"][idx])
    p4 = np.column_stack([data["gamma"][idx] * mass,
                          (data["beta_gamma"][idx] * mass)[:, None] * direction])
    d, passed, _ = build_event_mc(p4, direction, entry_d[idx], exit_d[idx],
                                  templates, decay_samples, rng)
    _, N = scan_u2(d, passed, exit_d[idx] - entry_d[idx], w,
                   data["beta_gamma"][idx], float(templates["ctau_m_u2eq1"]),
                   L_INT_PB, u2_grid)
    return N


def run(paths: ModelPaths, flavors, masses, delta_bc, max_hit_events, decay_samples,
        seed_salt=""):
    mesh = get_mesh()
    u2_grid = np.logspace(LOG_U2_MIN, LOG_U2_MAX, N_U2_POINTS)
    rows = []
    for flavor in flavors:
        for mass in masses:
            ml = format_mass_for_filename(mass)
            templates = load_template_bundle(paths.templates / flavor / f"templates_{ml}.npz")
            if templates is None:
                continue
            base = paths.vectors / flavor
            N_tot = _scan_csv(base / "combined" / f"mN_{ml}.csv", flavor, mass, "comb",
                              mesh, templates, u2_grid, max_hit_events,
                              decay_samples, seed_salt)
            N_bc = _scan_csv(base / "Bc" / f"mN_{ml}.csv", flavor, mass, "Bc",
                             mesh, templates, u2_grid, max_hit_events,
                             decay_samples, seed_salt)

            def u2min(N):
                r = find_exclusion_band(u2_grid, N, N_THRESHOLD)
                return r["u2_min"] if r["has_sensitivity"] else np.nan

            c = u2min(N_tot)
            hi = u2min(N_tot + delta_bc * N_bc)   # Bc up   -> lower u2_min (stronger)
            lo = u2min(N_tot - delta_bc * N_bc)   # Bc down -> higher u2_min (weaker)
            f_bc = np.nan
            if np.isfinite(c) and c > 0:
                iu = int(np.argmin(np.abs(u2_grid - c)))
                f_bc = N_bc[iu] / N_tot[iu] if N_tot[iu] > 0 else 0.0
            rows.append({"flavor": flavor, "mass_GeV": mass, "delta_bc": delta_bc,
                         "u2_min": c, "u2_min_bc_lo": lo, "u2_min_bc_hi": hi,
                         "bc_frac_at_u2min": f_bc,
                         "decay_samples": decay_samples,
                         "max_hit_events": max_hit_events,
                         "seed_salt": seed_salt})
            shift = (np.log10(lo / hi) if np.isfinite(lo) and np.isfinite(hi) and hi > 0
                     else np.nan)
            print(f"  {flavor:5s} m={mass:<4} Bc_frac={f_bc*100:4.0f}%  "
                  f"u2_min={c:.3e}  band=[{hi:.3e},{lo:.3e}]  width={shift:.3f} dex",
                  flush=True)
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flavor", nargs="+", default=["Ue", "Umu", "Utau"])
    ap.add_argument("--mass", nargs="+", type=float, default=None)
    ap.add_argument("--delta-bc", type=float, default=SIGMA_BC_REL_UNCERT)
    ap.add_argument("--max-hit-events", type=int, default=4000)
    ap.add_argument("--decay-samples", type=int, default=50)
    ap.add_argument("--seed-salt", default="",
                    help="deterministic salt for independent numerical-control repeats")
    ap.add_argument("--vectors-dir", type=Path, default=None)
    ap.add_argument("--templates-dir", type=Path, default=None)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    paths = ModelPaths.resolve("hnl", vectors=args.vectors_dir, templates=args.templates_dir)
    df = run(paths, args.flavor, args.mass or DEFAULT_MASSES, args.delta_bc,
             args.max_hit_events, args.decay_samples, args.seed_salt)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"wrote {args.out} ({len(df)} rows, delta_bc={args.delta_bc})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
