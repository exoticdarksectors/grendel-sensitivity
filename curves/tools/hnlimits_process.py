"""Write the transformed HNLimits regions the paper panels draw.

For each single-flavour scenario, every record of the HNLimits workbook
sheet is transformed with HNLimits' own rules (mass units, confidence
level, Dirac/Majorana) and written to
``curves/data/hnl/hnlimits_<scenario>_<flavor>_majorana/<id>_<flavor>_majorana.dat``
with its provenance in the header. The renderer computes the same regions
in memory (``curves.plot.hnl_panel.processed_curves``); these files are the
inspectable record of what it draws.

    python -m curves.tools.hnlimits_process            # rewrite the tracked files
    python -m curves.tools.hnlimits_process --check    # regenerate into a temporary
                                                       # directory and compare bytes
"""
from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

from ..plot import hnl_panel
from ..plot.data import HNL_DATA, HNLIMITS_DATABASE


def regenerate(root: Path) -> list[Path]:
    written = []
    for scenario, config in hnl_panel.SCENARIOS.items():
        records = hnl_panel.load_metadata(config["sheet"])
        out_dir = root / hnl_panel.processed_dir(scenario, config).name
        out_dir.mkdir(parents=True, exist_ok=True)
        for item in hnl_panel.processed_curves(scenario, config, records):
            path = out_dir / Path(item["path"]).name
            path.write_text(hnl_panel.processed_file_text(scenario, config, item))
            written.append(path)
    return written


def check(root: Path = HNL_DATA) -> int:
    with tempfile.TemporaryDirectory() as tmp:
        fresh = regenerate(Path(tmp))
        rel = {p.relative_to(tmp) for p in fresh}
        tracked = {p.relative_to(root) for d in root.glob("hnlimits_*_majorana") for p in d.glob("*.dat")}
        bad = 0
        for missing in sorted(tracked - rel):
            print(f"not regenerated: {missing}")
            bad += 1
        for extra in sorted(rel - tracked):
            print(f"not tracked:     {extra}")
            bad += 1
        for r in sorted(rel & tracked):
            if (Path(tmp) / r).read_bytes() != (root / r).read_bytes():
                print(f"differs:         {r}")
                bad += 1
    print(f"{len(rel & tracked)} files compared, {bad} problems (workbook {HNLIMITS_DATABASE.name})")
    return 1 if bad else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="compare a regeneration with the tracked files")
    args = ap.parse_args(argv)
    if args.check:
        return check()
    for path in regenerate(HNL_DATA):
        print(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
