"""BC5 Higgs production, stage 3: the tracked Higgs (pT, y) histogram and its sampler."""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

from .higgs_production import mg5_work_dir
from .model import M_HIGGS

POOL_PATH = Path(__file__).resolve().parent / "data" / "higgs_pt_y_14TeV.csv.gz"

SIGMA_H_YR4_PB = {
    "ggf": 54.61 + 0.5955,
    "vbf": 4.275,
    "wh": 1.510,
    "zh": 0.9836,
    "tth": 0.6128,
}

PT_EDGES = np.concatenate([np.arange(0.0, 10.0, 1.0), np.arange(10.0, 50.0, 2.5),
                           np.arange(50.0, 200.0, 10.0), np.arange(200.0, 500.0, 50.0),
                           np.arange(500.0, 1501.0, 250.0)])
Y_EDGES = np.round(np.arange(-6.0, 6.0 + 1e-9, 0.2), 6)


def build_pool(work: Path, out: Path = POOL_PATH, modes=None) -> Path:
    """Bin the showered samples of ``modes`` (default: every mode of ``SIGMA_H_YR4_PB``) into the
    histogram file ``out``."""
    shower = json.loads((work / "higgs_shower_manifest.json").read_text())
    lhe = json.loads((work / "higgs_lhe_manifest.json").read_text())
    modes = list(modes or SIGMA_H_YR4_PB)
    blocks, header = [], []
    header.append("Higgs (pT, y) histograms per production mode, pp at 14 TeV, for BC5 h -> S S")
    header.append(f"built {datetime.now().isoformat()} by grendel.models.scalar.higgs_pool")
    header.append("hard process: MadGraph5_aMC@NLO v3.6.6, LO, NNPDF40_nlo_as_01180 (LHAPDF 331700), cards/proc_card_<mode>.dat")
    header.append("shower: Pythia 8 ISR+FSR (no MPI, no hadronisation), Higgs stable; grendel.models.scalar.higgs_shower")
    header.append("normalisation used by the production: HXSWG YR4 14 TeV (m_H = 125.09 GeV), " +
                  ", ".join(f"{k}={v:.4g} pb" for k, v in SIGMA_H_YR4_PB.items()) +
                  "; ggf = ggF N3LO 54.61 + bbH 0.5955 (bbH carried with the ggF shape); tH omitted")
    header.append("columns: mode, pt_lo_GeV, pt_hi_GeV, y_lo, y_hi, count")
    for mode in modes:
        if mode not in shower:
            raise FileNotFoundError(f"no showered sample for mode {mode!r} in {work}")
        rows = np.loadtxt(shower[mode]["csv"], delimiter=",", comments="#").reshape(-1, 6)
        counts, _, _ = np.histogram2d(rows[:, 0], rows[:, 1], bins=[PT_EDGES, Y_EDGES])
        lost = len(rows) - counts.sum()
        header.append(f"{mode}: n_higgs={len(rows)} (outside grid: {int(lost)}), pythia={shower[mode]['pythia_version']}, "
                      f"seed={shower[mode]['seed']}, sigma_mg5_lo_pb={lhe[mode]['sigma_lo_pb']}, "
                      f"<pT>={rows[:, 0].mean():.2f} GeV, <|y|>={np.abs(rows[:, 1]).mean():.3f}")
        ipt, iy = np.nonzero(counts)
        for a, b in zip(ipt, iy):
            blocks.append(f"{mode},{PT_EDGES[a]:.6g},{PT_EDGES[a + 1]:.6g},{Y_EDGES[b]:.6g},{Y_EDGES[b + 1]:.6g},{int(counts[a, b])}")
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt") as fh:
        fh.write("".join(f"# {line}\n" for line in header))
        fh.write("\n".join(blocks) + "\n")
    return out


_POOL = None


def load_pool(path: Path = POOL_PATH):
    """{mode: (pt_lo, pt_hi, y_lo, y_hi, count)} arrays, cached."""
    global _POOL
    if _POOL is None or _POOL[0] != path:
        table = {}
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                mode, *values = line.strip().split(",")
                table.setdefault(mode, []).append([float(v) for v in values])
        _POOL = (path, {mode: np.array(rows, float) for mode, rows in table.items()})
    return _POOL[1]


def sample_higgs_4vectors(n, rng, *, sigma_pb=None, path: Path = POOL_PATH):
    """``n`` Higgs four-vectors drawn from the tracked histograms, the mode by its cross section, the
    cell by its count, uniform inside the cell, phi uniform."""
    sigma_pb = dict(sigma_pb or SIGMA_H_YR4_PB)
    pool = load_pool(path)
    modes = [m for m in sigma_pb if m in pool]
    weights = np.array([sigma_pb[m] for m in modes], float)
    total = float(weights.sum())
    mode_idx = rng.choice(len(modes), size=n, p=weights / total)
    pt = np.empty(n)
    y = np.empty(n)
    for i, mode in enumerate(modes):
        sel = np.flatnonzero(mode_idx == i)
        if len(sel) == 0:
            continue
        cells = pool[mode]
        prob = cells[:, 4] / cells[:, 4].sum()
        pick = rng.choice(len(cells), size=len(sel), p=prob)
        pt[sel] = rng.uniform(cells[pick, 0], cells[pick, 1])
        y[sel] = rng.uniform(cells[pick, 2], cells[pick, 3])
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    mt = np.sqrt(pt ** 2 + M_HIGGS ** 2)
    return {"E": mt * np.cosh(y), "px": pt * np.cos(phi), "py": pt * np.sin(phi),
            "pz": mt * np.sinh(y), "pt": pt, "y": y,
            "mode": np.array(modes, dtype=object)[mode_idx], "sigma_pb": total}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work-dir", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=POOL_PATH)
    ap.add_argument("--modes", nargs="+", default=None, help="default: every mode")
    args = ap.parse_args(argv)
    out = build_pool(args.work_dir or mg5_work_dir(), args.out, modes=args.modes)
    pool = load_pool(out)
    for mode, cells in pool.items():
        print(f"{mode}: {len(cells)} cells, {int(cells[:, 4].sum())} Higgs")
    print(f"Saved: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
