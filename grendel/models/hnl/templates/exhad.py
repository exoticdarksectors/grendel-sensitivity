"""exHad HNL rest-frame decay templates."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from ....io.paths import ModelPaths
from ....reco.exhad_generation import (TemplateModel, TemplatePoint, add_generation_options,
                                        default_exhad_root, run)
from ....reco.templates import HBARC_GEV_M, import_exhad
from ..mass_grid import MASS_GRID, format_mass_for_filename

FLAVOR_TO_MIXING = {"Ue": (1.0, 0.0, 0.0), "Umu": (0.0, 1.0, 0.0), "Utau": (0.0, 0.0, 1.0)}
WIDTH_COLUMN = {"Ue": 1, "Umu": 2, "Utau": 3}
SEED_NAMESPACE = "exhad-hnl-templates/v1"


def load_width_table(exhad_root: Path) -> np.ndarray:
    """exHad HNL total widths at |U_alpha|^2 = 1: columns m_N, Gamma_e, Gamma_mu, Gamma_tau [GeV]."""
    model_info = import_exhad(exhad_root).model_info
    path = Path(model_info("hnl")["tables"]["widths"])
    table = np.loadtxt(path)
    if table.ndim != 2 or table.shape[1] != 4:
        raise RuntimeError(f"unexpected width table shape {table.shape} in {path}")
    return table


def ctau_u2eq1_exhad(mass_GeV: float, flavor: str, table: np.ndarray) -> float:
    """Proper decay length at U^2 = 1 in metres (log-log interpolation)."""
    m = table[:, 0]
    g = table[:, WIDTH_COLUMN[flavor]]
    if not (m[0] <= mass_GeV <= m[-1]):
        raise ValueError(f"m_N = {mass_GeV} GeV outside the exHad width table [{m[0]}, {m[-1]}]")
    log_g = np.interp(np.log(mass_GeV), np.log(m), np.log(g))
    return HBARC_GEV_M / float(np.exp(log_g))


def _reference_ctau(reference_dir, pt: TemplatePoint) -> float:
    if reference_dir is None:
        return np.nan
    path = Path(reference_dir) / pt.flavor / f"templates_{pt.label}.npz"
    if not path.exists():
        return np.nan
    with np.load(path) as ref:
        return float(ref["ctau_m_u2eq1"])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flavor", nargs="+", default=["Ue", "Umu", "Utau"], choices=list(FLAVOR_TO_MIXING))
    ap.add_argument("--mass", type=float, nargs="+", default=None,
                    help="restrict to these masses (default: the full grid)")
    add_generation_options(ap, robust=False)
    ap.add_argument("--reference-templates", type=Path, default=None,
                    help="FairShip template dir; its ctau is stored as ctau_m_u2eq1_reference")
    args = ap.parse_args(argv)
    if args.n_templates < 1:
        ap.error("--n-templates must be >= 1")

    width_table = load_width_table(args.exhad_root.expanduser().resolve())
    model = TemplateModel(
        exhad_model="hnl", decay_model="exhad", flavor_tag=None,
        seed_namespace=SEED_NAMESPACE,
        ctau_source="exhad:data/hnl/eventcalc_total_widths.dat (arXiv:1805.08567), log-log interpolated",
        ctau_ref=lambda pt: ctau_u2eq1_exhad(pt.mass, pt.flavor, width_table),
        reference_ctau=lambda pt: _reference_ctau(args.reference_templates, pt),
        dest=lambda pt, out: out / pt.flavor / f"templates_{pt.label}.npz",
        strategy="attempt",
        manifest_extra={"flavors": args.flavor},
    )
    masses = args.mass if args.mass is not None else list(MASS_GRID)
    points = [TemplatePoint(float(m), format_mass_for_filename(m), fl, FLAVOR_TO_MIXING[fl])
              for fl in args.flavor for m in masses]
    out_dir = args.out or ModelPaths.resolve("hnl").templates
    return run(model, points, out_dir, args)


if __name__ == "__main__":
    sys.exit(main())
