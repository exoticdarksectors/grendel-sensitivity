"""Event-level check of the four-vector pools through the ray-cast caches of the scans.

The published pool CSVs are not readable here (validate_pools.py), but bc10/maksym_alp_fix_check_20261018 kept the
ray-cast caches that the code of the v2.3/v2.4 scans (c621279, UP542 tunnel 2026-10-16.1) computed on macOS from the
v2.2 pool for 18 masses: per four-vector, whether it hits the fiducial volume and its entry and exit distances.
Our scans write the same caches (04_scan.sh keeps them). Same hits arrays and entry/exit distances equal to rounding
mean the same events in the same order with the same directions.

    python validate_geometry.py <work/alp/scan_VARIANT> <report.txt>

The stored caches are reference material ($GRENDEL_REF); without them the script prints SKIP and exits 0.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from ref import REF as REF_TREE, skip_unless  # noqa: E402  (pipelines/common/ref.py)

REF = REF_TREE / "grendel/campaigns/bc10/maksym_alp_fix_check_20261018/work/grendel_bc10/geometry_cache"


def main(scan_dir, report):
    skip_unless("BC10 ray-cast caches vs the stored macOS ones", REF)
    scan_dir = Path(scan_dir)
    lines, n_ok, n_cmp = [], 0, 0
    for ref in sorted(REF.glob("geom_*.npz")):
        lab = ref.stem[5:]
        mass = float(lab.replace("p", "."))
        # scan directories are named as in points_all.txt: "1.0", "1.5", "3.42" (not "1")
        names = dict.fromkeys([f"{mass:g}", f"{mass:.1f}"])
        ours = next((p for n in names for p in sorted(scan_dir.glob(f"bc10_{n}/geo/{ref.name}"))), None)
        if ours is None:
            lines.append(f"{mass:g}: no cache of ours yet")
            continue
        a, b = np.load(ours), np.load(ref)
        n_cmp += 1
        ha, hb = a["hits"].astype(bool), b["hits"].astype(bool)
        if ha.shape != hb.shape:
            lines.append(f"{mass:g}: {ha.size} vs {hb.size} four-vectors")
            continue
        both = ha & hb
        de = np.nanmax(np.abs(a["entry_d"][both] - b["entry_d"][both])) if both.any() else 0.0
        dx = np.nanmax(np.abs(a["exit_d"][both] - b["exit_d"][both])) if both.any() else 0.0
        flips = int((ha != hb).sum())
        ok = flips == 0 and max(de, dx) < 1e-6
        n_ok += ok
        lines.append(f"{mass:g}: {ha.size} four-vectors, hits {int(ha.sum())} / stored {int(hb.sum())}, "
                     f"hit flags differing {flips}, max |d entry| {de:.1e} m, max |d exit| {dx:.1e} m"
                     + ("" if ok else "  DIFFER"))
    head = f"ray-cast caches compared with the stored macOS ones: {n_ok} of {n_cmp} identical events (of {len(list(REF.glob('geom_*.npz')))} stored)"
    Path(report).write_text(head + "\n" + "\n".join(lines) + "\n")
    print(head)
    return 0 if n_ok == n_cmp else 1


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
