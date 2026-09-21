"""exHad BC4 dark-scalar rest-frame decay templates.

Replaces the analytic two-body decay engine of ``acceptance`` (every mode as
``S -> x xbar`` back to back; ``s s``, ``c c``, ``g g`` and ``4 pi`` as two
leading charged pions) with exHad's full final states (Kryshtal &
Ovchynnikov, arXiv:2609.16104): separately generated ``pi pi``, ``K K``,
``p pbar``, ``n nbar`` and four-pion channels below 2 GeV, and above 2 GeV
the constrained string fragmentation of the gluon-pair and ``s sbar``
widths. The default model ``scalar-1809`` takes the widths of
arXiv:1809.01876 with the corrections and channels of arXiv:1904.10447;
``scalar-central``/``-lower``/``-upper`` use the dispersive calculation of
arXiv:2407.13587.

Lifetime: ``ctau_m_u2eq1`` is exHad's own table for the chosen model, so
the channel mix and the decay probability come from one width model. It
differs from the Winkler digitisation behind the first published curve by
up to a factor 2 above 2 GeV (gluon-pair width); the reference is stored
alongside so the two effects can be separated.

    python -m grendel.models.scalar.templates_exhad --n-templates 20000
"""
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
