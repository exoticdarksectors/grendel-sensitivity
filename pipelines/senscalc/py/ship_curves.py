"""SHiP BC10 curve files.

    python -I ship_curves.py <out dir> <Maksym JSON | -> <rerun JSON>

1. Regenerates the tracked curves/data/bc10/projections/ship_senscalc.dat from M. Ovchynnikov's JSON with the
   repository tool itself (curves.tools.senscalc_bc10_ship, unchanged since v2.4.0), redirected to
   <out dir>/projections/ship_senscalc.dat (the source JSON is linked into a private GRENDEL_CURVES_SOURCES).
   The JSON is a private communication, kept in the optional reference tree only; "-" skips this part.
2. Converts our headless SHiP-ECN3 rerun (same pipeline as the FASER2 curve) to the same format,
   <out dir>/projections/ship_senscalc_rerun.dat. This is a validation product, not a replacement: the paper
   draws the author's curve.
Must run from src/grendel-sensitivity with the analysis venv (imports curves.*).
"""
import json
import os
import sys
from pathlib import Path


def main() -> int:
    out, maksym, rerun = Path(sys.argv[1]).resolve(), sys.argv[2], Path(sys.argv[3])
    src = out / "curves_sources"
    rel = Path("arXiv-2501.04525/Sensitivity_ALP-fermion_at_SHiP-ECN3_Nev=2.3_Npot=6.e20.json")
    (src / rel).parent.mkdir(parents=True, exist_ok=True)
    if (src / rel).is_symlink() or (src / rel).exists():
        (src / rel).unlink()
    if maksym != "-":
        (src / rel).symlink_to(Path(maksym))
    os.environ["GRENDEL_CURVES_SOURCES"] = str(src)

    import curves.tools.senscalc_bc10_ship as tool  # noqa: E402  (reads GRENDEL_CURVES_SOURCES at import)
    from curves.tools import pdfio  # noqa: E402

    if maksym != "-":
        tool.BC10_DATA = out  # writes out/projections/ship_senscalc.dat
        tool.main()

    variants = json.loads(rerun.read_text())
    central = [v for v in variants if all(v[k] == val for k, val in tool.CENTRAL.items())]
    assert len(central) == 1, len(central)
    pdfio.write_chains(
        out / "projections" / "ship_senscalc_rerun.dat",
        [[(m, gy / tool.V_H) for m, gy in dom] for dom in central[0]["Sensitivity domains"]],
        ["SHiP fermiophilic-ALP (BC10) sensitivity, SHiP-ECN3 rerun with SensCalc eb2c456 headless on Linux",
         "(grendel-repro pipelines/senscalc, the pipeline of campaigns/bc10/faser2_senscalc_20261025), central",
         "bremsstrahlung and Drell-Yan variant; settings of ship_senscalc.dat (6e20 PoT, N_ev >= 2.3).",
         f"Source file: {rerun.name}",
         "Source y axis g_Y = 2 v/f_a; stored as 1/f (BNT) = g_Y/246. Validation product, not drawn."],
        columns='Format: chains separated by blank line; each row "mass_GeV  invf_GeV^-1"')
    return 0


if __name__ == "__main__":
    sys.exit(main())
