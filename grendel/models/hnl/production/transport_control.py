#!/usr/bin/env python3
"""Transport control: is the displaced charged-kaon decay origin immaterial?

The default kaon model (`generate_kaon_csvs.py`) treats charged-kaon transport as
a per-kaon survival weight and casts the HNL from the IP, on the stated grounds
that the displaced decay origin (<= `d_esc` ~ 1.5 m) does not change the HNL
detector acceptance. This script is the tracked evidence for that claim. It
quantifies the approximation against the real GRENDEL fiducial mesh:

  1. sample charged kaons from the committed Pythia SoftQCD spectrum;
  2. for each kaon sample a decay point along its flight (truncated-exponential
     within `d_esc`) and carry the transport survival weight `w`;
  3. decay K -> l N to an HNL using the production kinematics
     (`_sample_hnl_from_mesons`);
  4. for a scan of HNL proper lifetimes, ray-cast each HNL against the fiducial
     mesh from BOTH the IP and its displaced kaon-decay point, weight by the
     decay-in-volume probability `exp(-d_in/L) - exp(-d_out/L)` (L = beta*gamma*ctau,
     the same form the sensitivity engine uses), and compare the summed accepted
     yields.

The reported ratio (displaced / IP) per proper lifetime is the number quoted as
"~1.0" in `REMAINING_WORK.md` and `production/constants.py`. Results are written
to `production/data/transport_control.json`.

    python -m production.decay_engine.transport_control            # default scan
    python production/decay_engine/transport_control.py --n 300000
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

from grendel.production.constants import M_KAON, LEPTON_MASSES, FLAVOR_TO_LEPTON_PDG, KAON_D_ESC
from grendel.models.hnl.production.generate_kaon_csvs import (
    sample_kaon_4vectors_pythia, CTAU_KAON, KAON_PDG,
)
from grendel.models.hnl.production.generate_meson_csvs import (
    compute_production_br_components, _sample_hnl_from_mesons,
)
from grendel.models.hnl.production.hnlcalc import init_hnlcalc

DATA = Path(__file__).resolve().parent / "data" / "transport_control.json"
# HNL proper lifetimes (m) spanning short-lived (decays before the y~22 m detector)
# through the long-lifetime plateau that sets the sensitivity edge.
CTAU_SCAN_M = [0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0]


def _unit(px, py, pz):
    p = np.sqrt(px * px + py * py + pz * pz)
    return np.stack([px / p, py / p, pz / p], axis=1), p


def cast(mesh, origins, directions):
    """Ray-cast per-ray (origin, direction) against the fiducial mesh, returning
    hit mask and entry/exit distances measured FROM each ray's own origin. Mirrors
    analysis._engine.compute_geometry but allows a distinct origin per ray."""
    n = len(directions)
    hits = np.zeros(n, dtype=bool)
    entry_d = np.full(n, np.nan)
    exit_d = np.full(n, np.nan)
    cand = np.where(directions[:, 1] > 0.01)[0]  # +y half-space (same as the engine)
    CHUNK = 25000
    for cs in range(0, len(cand), CHUNK):
        ci = cand[cs:cs + CHUNK]
        loc, rid, _ = mesh.ray.intersects_location(
            ray_origins=origins[ci], ray_directions=directions[ci])
        if not len(loc):
            continue
        d = np.linalg.norm(loc - origins[ci][rid], axis=1)
        order = np.argsort(rid)
        rid_s, d_s = rid[order], d[order]
        uniq, start, cnt = np.unique(rid_s, return_index=True, return_counts=True)
        for j in np.where(cnt >= 2)[0]:
            gi = ci[uniq[j]]
            seg = np.sort(d_s[start[j]:start[j] + cnt[j]])
            hits[gi] = True
            entry_d[gi], exit_d[gi] = seg[0], seg[1]
    return hits, entry_d, exit_d


def _decay_in_volume(hits, entry_d, exit_d, lam):
    p = np.zeros(len(hits))
    m = hits & np.isfinite(entry_d) & np.isfinite(exit_d)
    p[m] = np.exp(-entry_d[m] / lam[m]) - np.exp(-exit_d[m] / lam[m])
    return p


def run(n_kaons, m_N, flavor, d_esc, seed):
    rng = np.random.default_rng(seed)
    pool = sample_kaon_4vectors_pythia(n_kaons, rng)
    kdir, p_k = _unit(pool["px"], pool["py"], pool["pz"])
    lam_k = (p_k / M_KAON) * CTAU_KAON                     # kaon lab decay length
    w = 1.0 - np.exp(-d_esc / lam_k)                       # transport survival weight
    # decay point along the flight, conditioned on decaying within d_esc
    u = rng.random(n_kaons)
    length = -lam_k * np.log(1.0 - u * (1.0 - np.exp(-d_esc / lam_k)))
    displaced = kdir * length[:, None]                    # m, from the IP

    hnl = init_hnlcalc(flavor)
    lepton_pdg = FLAVOR_TO_LEPTON_PDG[flavor]
    m_lepton = LEPTON_MASSES[flavor]
    _, br_3body, br_total = compute_production_br_components(hnl, KAON_PDG, lepton_pdg, m_N)
    hnl4 = _sample_hnl_from_mesons(
        pool["E"], pool["px"], pool["py"], pool["pz"],
        M_KAON, m_lepton, m_N, br_3body, br_total, hnl, rng)
    ndir, p_N = _unit(hnl4[:, 1], hnl4[:, 2], hnl4[:, 3])
    bg_N = p_N / m_N

    from grendel.constants import CMS_ORIGIN
    from grendel.geometry.raycast import get_mesh
    mesh = get_mesh()
    ip = np.tile(np.asarray(CMS_ORIGIN, dtype=float), (n_kaons, 1))
    hits_ip, ein_ip, eout_ip = cast(mesh, ip, ndir)
    hits_dp, ein_dp, eout_dp = cast(mesh, displaced, ndir)

    rows = []
    for ctau in CTAU_SCAN_M:
        lam = bg_N * ctau
        P_ip = _decay_in_volume(hits_ip, ein_ip, eout_ip, lam)
        P_dp = _decay_in_volume(hits_dp, ein_dp, eout_dp, lam)
        den = float(np.sum(w * P_ip))
        num = float(np.sum(w * P_dp))
        rows.append({
            "ctau_m": ctau,
            "yield_ip": den,
            "yield_displaced": num,
            "ratio_displaced_over_ip": (num / den) if den > 0 else float("nan"),
        })
    return {
        "config": {
            "n_kaons": n_kaons, "m_N_GeV": m_N, "flavor": flavor,
            "d_esc_m": d_esc, "ctau_kaon_m": CTAU_KAON, "seed": seed,
            "mean_transport_weight": float(np.mean(w)),
            "median_displacement_m": float(np.median(length)),
            "n_hit_ip": int(hits_ip.sum()), "n_hit_displaced": int(hits_dp.sum()),
        },
        "scan": rows,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n", type=int, default=300_000, help="number of kaons to sample")
    ap.add_argument("--m-n", type=float, default=0.2, help="HNL mass (GeV), in the K->lN window")
    ap.add_argument("--flavor", default="Umu", choices=["Ue", "Umu"])
    ap.add_argument("--d-esc", type=float, default=KAON_D_ESC)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("-o", "--output", type=Path, default=DATA)
    args = ap.parse_args(argv)

    result = run(args.n, args.m_n, args.flavor, args.d_esc, args.seed)
    c = result["config"]
    print(f"m_N={c['m_N_GeV']} GeV {c['flavor']}, {c['n_kaons']} kaons, "
          f"d_esc={c['d_esc_m']} m, <w>={c['mean_transport_weight']:.3f}, "
          f"median displacement={c['median_displacement_m']:.3f} m", file=sys.stderr)
    print(f"hits: IP={c['n_hit_ip']}  displaced={c['n_hit_displaced']}", file=sys.stderr)
    print(f"{'ctau_N [m]':>12} {'ratio disp/IP':>16}", file=sys.stderr)
    for r in result["scan"]:
        print(f"{r['ctau_m']:>12.1f} {r['ratio_displaced_over_ip']:>16.4f}", file=sys.stderr)

    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
