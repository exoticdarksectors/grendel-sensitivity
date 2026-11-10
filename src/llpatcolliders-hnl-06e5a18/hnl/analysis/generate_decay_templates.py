"""Generate FairShip HNL rest-frame decay templates and cache them to .npz.

This is the ROOT+Pythia8 stage of the signal-acceptance pipeline (see
``hnl/README.md`` and ``decay_reco_acceptance.py``). It must run under a Python
that has PyROOT + Pythia8 (``ROOT.TPythia8``); on this machine that is the
Homebrew-ROOT venv, NOT the conda ``hnl`` env. The downstream analysis consumes
the cached templates in the conda env without ROOT.

For each (flavor, mass) it samples ``--n-templates`` rest-frame decays with
FairShip and writes a ragged-but-flat array bundle so the analysis can split
per template and boost to the lab:

    out/<flavor>/templates_<masslabel>.npz
        daughter_counts (N,)            int32   daughters per template
        pdg             (M,)            int32   M = sum(daughter_counts)
        px, py, pz, energy, mass (M,)   float64 rest-frame 4-momentum + mass
        charge          (M,)            float64
        stable          (M,)            bool    final-state (status == 1)
        mass_GeV, ctau_m_u2eq1, n_templates, seed   scalars
        flavor                          str

All daughters (including invisible decays) are stored, so the visible
branching fraction is recoverable as the fraction of templates with >= 2
charged stable daughters.

Example (run with the FairShip venv python):

    /Volumes/sandbox/projects/aaaPHYSICSaaa/.venvs/fairship/bin/python \
        hnl/analysis/generate_decay_templates.py \
        --flavor Ue Umu Utau --n-templates 20000 --out hnl/tmp/decay_templates
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_ANALYSIS_DIR = Path(__file__).resolve().parent
_HNL_ROOT = _ANALYSIS_DIR.parent
for _p in (str(_ANALYSIS_DIR), str(_HNL_ROOT)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import fairship_decay as fd  # noqa: E402  (needs PyROOT + Pythia8)
from config_mass_grid import MASS_GRID, format_mass_for_filename  # noqa: E402


def _flatten_templates(templates):
    """Pack a list of DecayTemplate into flat arrays + per-template counts."""
    counts = np.array([len(t.pdg) for t in templates], dtype=np.int32)
    if counts.sum() == 0:
        empty_f = np.empty(0, dtype=np.float64)
        return dict(
            daughter_counts=counts,
            pdg=np.empty(0, dtype=np.int32),
            px=empty_f, py=empty_f, pz=empty_f, energy=empty_f, mass=empty_f,
            charge=empty_f, stable=np.empty(0, dtype=bool),
        )
    cat = lambda attr, dt: np.concatenate(
        [np.asarray(getattr(t, attr), dtype=dt) for t in templates])
    return dict(
        daughter_counts=counts,
        pdg=cat("pdg", np.int32),
        px=cat("px", np.float64), py=cat("py", np.float64),
        pz=cat("pz", np.float64), energy=cat("energy", np.float64),
        mass=cat("mass", np.float64), charge=cat("charge", np.float64),
        stable=cat("stable", bool),
    )


def generate_one(mass_GeV, flavor, n_templates, out_dir, seed, skip_existing):
    label = format_mass_for_filename(mass_GeV)
    dest = out_dir / flavor / f"templates_{label}.npz"
    if skip_existing and dest.exists():
        print(f"  [skip] {flavor} m={mass_GeV} ({dest.name})")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    templates = fd.sample_rest_frame_decays(
        mass_GeV, n_templates, flavor=flavor, random_seed=seed)
    bundle = _flatten_templates(list(templates))
    ctau_m = float(fd.ctau_u2eq1(mass_GeV, flavor=flavor))
    n_charged = np.array(
        [int(((np.abs(np.asarray(t.charge)) > 0.5) & np.asarray(t.stable)).sum())
         for t in templates])
    vis_frac = float((n_charged >= 2).mean()) if len(n_charged) else 0.0
    np.savez_compressed(
        dest,
        mass_GeV=np.float64(mass_GeV),
        ctau_m_u2eq1=np.float64(ctau_m),
        n_templates=np.int32(n_templates),
        seed=np.int64(-1 if seed is None else seed),
        flavor=np.array(flavor),
        **bundle,
    )
    print(f"  [ok]   {flavor} m={mass_GeV:.3f}  vis(>=2 chg)={vis_frac:5.1%}  "
          f"ctau(U2=1)={ctau_m:.3e} m  -> {dest.name}")
    return dest


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flavor", nargs="+", default=["Ue", "Umu", "Utau"],
                    choices=["Ue", "Umu", "Utau"])
    ap.add_argument("--mass", type=float, nargs="+", default=None,
                    help="restrict to these masses (default: full MASS_GRID)")
    ap.add_argument("--n-templates", type=int, default=20_000)
    ap.add_argument("--out", type=Path, default=_HNL_ROOT / "tmp" / "decay_templates")
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--skip-existing", action="store_true",
                    help="do not regenerate templates whose .npz already exists")
    args = ap.parse_args(argv)

    masses = args.mass if args.mass is not None else MASS_GRID
    out_dir = args.out.expanduser().resolve()
    print(f"Generating decay templates -> {out_dir}")
    print(f"  flavors={args.flavor}  masses={len(masses)}  "
          f"n_templates={args.n_templates}  seed={args.seed}")
    for flavor in args.flavor:
        for m in masses:
            generate_one(m, flavor, args.n_templates, out_dir, args.seed,
                         args.skip_existing)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
