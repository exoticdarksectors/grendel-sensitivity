"""Raw HNLimits curves for the three single-flavour scenarios.

Copies the per-experiment existing-exclusion files listed in each
scenario sheet of the HNLimits workbook (``Ue4``, ``Umu4``, ``Utau4``) into
``curves/data/hnl/hnlimits_{Ue,Umu,Utau}/``, and the projection files into
``hnlimits_projections_{flavor}/``. HNLimits exposes only the Umu
projections through its sheet; the Ue and Utau ones exist as bare ``.dat``
files inside ``include/data/<arxiv_id>/`` and are listed here by hand.

These are the untransformed curves (no unit, CL or Dirac/Majorana
conversion); the diagnostic figures read them. The transformed regions the
paper panels draw are made by ``hnlimits_process``.

    python -m curves.tools.hnlimits_pull
"""
from __future__ import annotations

import sys

import numpy as np

from ..plot import hnl_panel
from ..plot.data import HNL_DATA

PROJECTIONS = {
    "Ue": [("SHiP", "SHiP", "Ue4.dat", "Ue4_top.dat", "arXiv:1811.00930", "Current SHiP"),
           ("FASER2", "FASER_2018", "converted_HNLe_contour_FASER2.txt", None, "arXiv:1811.12522",
            "FASER2 closed polygon")],
    "Umu": [("SHiP", "SHiP", "Umu4.dat", "Umu4_top.dat", "arXiv:1811.00930", "Current SHiP"),
            ("SHiP_TP_optimistic", "1811.00930", "Umu4_optimistic.dat", None, "arXiv:1811.00930",
             "SHiP TP optimistic"),
            ("SHiP_TP_conservative", "1811.00930", "Umu4_conservative.dat", None, "arXiv:1811.00930",
             "SHiP TP conservative"),
            ("MATHUSLA", "1805.08567", "mathusla_lower.dat", "mathusla_upper.dat", "arXiv:1805.08567",
             "Bondarenko digitization of MATHUSLA paper"),
            ("FASER2", "FASER_2018", "converted_HNLmu_contour_FASER2.txt", None, "arXiv:1811.12522",
             "FASER2 closed polygon"),
            ("DUNE", "DUNE_Ballett", "Umu4.dat", "Umu4_top.dat", "arXiv:1905.00284",
             "DUNE near detector (Ballett et al.)"),
            ("FCCee", "1805.08567", "FCCee_lower.dat", "FCCee_upper.dat", "arXiv:1805.08567", "FCC-ee bonus")],
    "Utau": [("SHiP", "SHiP", "Utau4.dat", "Utau4_top.dat", "arXiv:1811.00930", "Current SHiP"),
             ("FASER2", "FASER_2018", "converted_HNLtau_contour_FASER2.txt", None, "arXiv:1811.12522",
              "FASER2 closed polygon")],
}


def load(rel):
    p = hnl_panel.hnlimits_data_root() / rel
    if not p.exists():
        return None
    a = np.loadtxt(p, comments="#")
    if a.ndim == 1:
        a = a.reshape(-1, 2)
    return a


def write_two_section(path, header, bottom, top, bottom_tag, top_tag):
    with path.open("w") as fh:
        for line in header:
            fh.write(f"# {line}\n")
        fh.write("# mass_GeV  U2\n")
        if bottom is not None:
            fh.write(f"# {bottom_tag}\n")
            for m, u2 in bottom:
                fh.write(f"{m:.5e}  {u2:.5e}\n")
        if top is not None:
            fh.write(f"\n# {top_tag}\n")
            for m, u2 in top:
                fh.write(f"{m:.5e}  {u2:.5e}\n")


def pull_projections():
    for flavor, items in PROJECTIONS.items():
        out_dir = HNL_DATA / f"hnlimits_projections_{flavor}"
        out_dir.mkdir(parents=True, exist_ok=True)
        for name, sub, bot_f, top_f, cite, note in items:
            b = load(f"{sub}/{bot_f}")
            t = load(f"{sub}/{top_f}") if top_f else None
            if b is None and t is None:
                print(f"  {flavor}/{name}: MISSING")
                continue
            path = out_dir / f"{name}_{flavor}.dat"
            write_two_section(path, [f"{name} HNL |U{flavor[1:] or 'e'}|^2 projection",
                                     f"Flavor scenario: single-flavor dominance in {flavor}",
                                     f"Source: {cite}  -- {note}",
                                     f"Pulled from HNLimits bundle ({sub}/{bot_f})",
                                     "Format: blank line separates lower edge from upper edge."],
                              b, t, "lower edge", "upper edge")
            print(f"  {flavor}/{name}: wrote {path}")


def pull_exclusions():
    for scenario, config in hnl_panel.SCENARIOS.items():
        flavor, sheet = config["flavor"], config["sheet"]
        out_dir = HNL_DATA / f"hnlimits_{flavor}"
        out_dir.mkdir(parents=True, exist_ok=True)
        n = 0
        for record in hnl_panel.load_metadata(sheet):
            b = load(record["file_bottom"]) if isinstance(record.get("file_bottom"), str) else None
            t = load(record["file_top"]) if isinstance(record.get("file_top"), str) else None
            if b is None and t is None:
                continue
            write_two_section(out_dir / f"{record['id']}_{flavor}.dat",
                              [f"{record['plot_label']}  (HNLimits id: {record['id']})",
                               f"Source: mhostert/Heavy-Neutrino-Limits sheet {sheet}",
                               f"Flavor: {flavor}",
                               'Format: blank line separates "bottom" and "top" curves.'],
                              b, t, "bottom", "top")
            n += 1
        print(f"  {flavor} exclusions: wrote {n} files to {out_dir}")


def main() -> int:
    pull_projections()
    pull_exclusions()
    return 0


if __name__ == "__main__":
    sys.exit(main())
