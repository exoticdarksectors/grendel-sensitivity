"""exHad BC10 fermiophilic-ALP rest-frame decay templates.

Full final states from exHad's ``alp-fermion`` model (Kryshtal &
Ovchynnikov, arXiv:2609.16104), replacing the Pythia-decayed exclusive
channels of ``templates_pythia``. exHad tabulates ctau at its own reference
coupling; ``Gamma(a -> mu mu)`` is convention-only, so the ratio of exHad's
to this package's mu mu width at masses where no hadronic channel is open
converts the table to the scan's reference coupling (the factor must be the
same at all four check masses to 1e-3).

    python -m grendel.models.alp.templates_exhad --n-templates 20000
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from ...io.paths import ModelPaths
from ...io.tables import open_table
from ...reco.exhad_generation import TemplateModel, TemplatePoint, add_generation_options, run
from ...reco.templates import HBARC_GEV_M, import_exhad, loglog_interp
from ...io.vectors import format_mass_for_filename
from . import model as alp_model
from .exclusive_decays import INV_F_REF, ctau_at_reference_coupling
from .mass_grid import ALP_MASS_GRID

SEED_NAMESPACE = "exhad-bc10-templates/v1"
CONVENTION_MASSES = (0.25, 0.30, 0.40, 0.50)   # only e e / mu mu / gamma gamma open


def load_exhad_tables(exhad_root: Path):
    """(ctau table (m, ctau) at exHad's reference coupling, mu mu branching-ratio nodes)."""
    info = import_exhad(exhad_root).model_info("alp-fermion")
    ctau = np.loadtxt(info["tables"]["ctau"])
    decay = json.loads(Path(info["tables"]["decay"]).read_text())
    mumu = next(entry for entry in decay if entry[0] == "muPmuM")
    return ctau, np.asarray(mumu[2], dtype=float)


def grendel_mumu_width(m_a: float) -> float:
    """Gamma(a -> mu mu) at INV_F_REF from the SensCalc export the scan uses."""
    data_dir = alp_model.SENSCALC_2501_DATA_DIR
    meta = json.loads((data_dir / "widths_metadata.json").read_text())
    column = next(c for c in meta["columns"] if c["canonical_name"] == "mumu")["source_index"] - 1
    table = np.loadtxt(open_table(data_dir / "widths_bnt.csv"), delimiter=",", skiprows=1)
    return float(np.interp(m_a, table[:, 0], table[:, column])) * INV_F_REF ** 2


def convention_factor(ctau_table, mumu_nodes) -> float:
    """ctau multiplier that brings exHad's table to the scan's reference coupling."""
    factors = []
    for m in CONVENTION_MASSES:
        ctau_ex = loglog_interp(m, ctau_table[:, 0], ctau_table[:, 1])
        br_mumu = float(np.interp(m, mumu_nodes[:, 0], mumu_nodes[:, 1]))
        gamma_ex = HBARC_GEV_M / ctau_ex * br_mumu
        factors.append(gamma_ex / grendel_mumu_width(m))
    factors = np.asarray(factors)
    spread = float(factors.max() / factors.min() - 1.0)
    if spread > 1e-3:
        raise RuntimeError(f"mu mu convention factor is not constant: {factors} (spread {spread:.2e})")
    return float(factors.mean())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mass", type=float, nargs="+", default=None, help="default: the BC10 mass grid")
    add_generation_options(ap, robust=True)
    args = ap.parse_args(argv)
    if args.n_templates < 1:
        ap.error("--n-templates must be >= 1")

    ctau_table, mumu_nodes = load_exhad_tables(args.exhad_root.expanduser().resolve())
    factor = convention_factor(ctau_table, mumu_nodes)
    print(f"  ctau convention factor exHad table -> INV_F_REF={INV_F_REF:g} GeV^-1: {factor:.4f}", flush=True)

    def extra_bundle(pt, bundle):
        commit = read_commit()
        return dict(
            includes_full_branching=np.bool_(True),
            primary_channel=bundle["channel_label"].astype("U32"),
            decay_backend=np.array(f"exHad-{commit[:12]}-alp-fermion-{args.variation}"),
            gluon_surrogate=np.array("none (exHad constrained fragmentation of the gg and ss widths)"),
            ctau_convention_factor=np.float64(factor),
            ctau_m_u2eq1_exhad_table=np.float64(loglog_interp(pt.mass, ctau_table[:, 0], ctau_table[:, 1])),
        )

    def read_commit():
        from ...reco.templates import read_exhad_provenance
        return read_exhad_provenance(args.exhad_root.expanduser().resolve())["exhad_commit"]

    masses = args.mass if args.mass is not None else list(ALP_MASS_GRID)
    tm = TemplateModel(
        exhad_model="alp-fermion", decay_model="exhad_alp_fermion_2501", flavor_tag="BC10",
        seed_namespace=SEED_NAMESPACE,
        ctau_source="exhad alp-fermion ctau table (log-log interpolated) x ctau_convention_factor",
        ctau_ref=lambda pt: loglog_interp(pt.mass, ctau_table[:, 0], ctau_table[:, 1]) * factor,
        reference_ctau=lambda pt: float(ctau_at_reference_coupling(pt.mass, "2501")),
        reference_source="exclusive_decays.ctau_at_reference_coupling(m, '2501'): SensCalc export behind the published scan",
        dest=lambda pt, out: out / f"templates_{pt.label}.npz",
        strategy="sub_request",
        extra_bundle=extra_bundle,
        skip=lambda pt: (lambda r: f"unsupported light-meson resonance ({r})" if r else None)(
            alp_model.excluded_light_meson_resonance(pt.mass)),
        manifest_extra={"ctau_convention_factor": factor,
                        "convention_masses_GeV": list(CONVENTION_MASSES),
                        "skipped_resonance_windows": [m for m in masses
                                                      if alp_model.excluded_light_meson_resonance(m)]},
    )
    points = [TemplatePoint(float(m), format_mass_for_filename(m)) for m in masses]
    out_dir = args.out or ModelPaths.resolve("bc10").templates
    return run(tm, points, out_dir, args)


if __name__ == "__main__":
    sys.exit(main())
