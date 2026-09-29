"""exHad BC4 dark-scalar rest-frame decay templates."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from ...io.paths import ModelPaths
from ...reco.exhad_generation import TemplateModel, TemplatePoint, add_generation_options, run
from ...reco.templates import import_exhad, loglog_interp
from . import model as smodel
from . import production as sprod

SCALAR_MODELS = ("scalar-1809", "scalar-central", "scalar-lower", "scalar-upper")
SEED_NAMESPACE = "exhad-bc4-templates/v1"


def load_ctau_table(exhad_root: Path, model_name: str) -> np.ndarray:
    """exHad (m_S [GeV], c*tau [m] at sin^2 theta = 1) table of a scalar model."""
    model_info = import_exhad(exhad_root).model_info
    path = Path(model_info(model_name)["tables"]["ctau"])
    table = np.asarray(json.loads(path.read_text()), dtype=float)
    if table.ndim != 2 or table.shape[1] != 2:
        raise RuntimeError(f"unexpected ctau table shape {table.shape} in {path}")
    return table


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="scalar-1809", choices=SCALAR_MODELS)
    ap.add_argument("--mass", type=float, nargs="+", default=None,
                    help="default: the BC4 mass grid")
    add_generation_options(ap, robust=True)
    args = ap.parse_args(argv)
    if args.n_templates < 1:
        ap.error("--n-templates must be >= 1")

    ctau_table = load_ctau_table(args.exhad_root.expanduser().resolve(), args.model)
    tm = TemplateModel(
        exhad_model=args.model, decay_model=f"exhad:{args.model}", flavor_tag="BC4",
        seed_namespace=SEED_NAMESPACE,
        ctau_source=f"exhad:{args.model} ctau table at sin^2 theta = 1, log-log interpolated",
        ctau_ref=lambda pt: loglog_interp(pt.mass, ctau_table[:, 0], ctau_table[:, 1]),
        reference_ctau=lambda pt: float(smodel.ctau_sin2theta1(pt.mass)),
        reference_source="scalar model ctau_sin2theta1 (Winkler arXiv:1809.01876 digitisation, first published curve)",
        dest=lambda pt, out: out / f"templates_{pt.label}.npz",
        strategy="sub_request",
    )
    masses = args.mass if args.mass is not None else list(sprod.MASS_GRID)
    points = [TemplatePoint(float(m), sprod._mass_label(m)) for m in masses]
    out_dir = args.out or ModelPaths.resolve("bc4").templates
    return run(tm, points, out_dir, args)


if __name__ == "__main__":
    sys.exit(main())
