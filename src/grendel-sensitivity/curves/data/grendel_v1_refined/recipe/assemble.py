"""Assemble the vector and template trees the two scans read.

v2: the new masses as produced (SIGMA_BC_PB = 1.64e6 pb) and the published v2
inputs at the check masses. v1: the new masses with the B_c rows rescaled to the
v1 rate (0.9e6 / 1.64e6, the inverse of review_fixes_20260928/make_hnl_vectors.py),
FairShip templates, and the published v1 inputs at the check masses.
"""
import os
import sys
from pathlib import Path

from grendel.io.vectors import format_mass_for_filename
from grendel.models.hnl.production.combine_channels import CHANNELS
from grendel.models.hnl.production.paths import llp_csv_path, read_csv_matrix, write_csv_matrix

import numpy as np

W = Path(os.environ.get("REFINE_DIR", Path(__file__).resolve().parent))   # the campaign directory
CAMP = W.parents[1]
NEW = [float(m) for m in sys.argv[1].split()]
CHECK = [float(m) for m in sys.argv[2].split()]
TAGS = sys.argv[3].split() if len(sys.argv) > 3 else ["v1", "v2"]
FLAVORS = ("Ue", "Umu", "Utau")
K_V1 = 0.9e6 / 1.64e6

PUBLISHED = {
    "v1": {"vectors": CAMP / "hnl/hnl_clean_20260825/llp_4vectors",
           "templates": CAMP / "hnl/hnl_exhad_20260917/tmp_fairship300/decay_templates"},
    "v2": {"vectors": CAMP / "review_fixes_20260928/hnl/llp_4vectors",
           "templates": CAMP / "v2_baseline_20260929/templates/hnl"},
}
NEW_TEMPLATES = {"v1": W / "templates_fairship", "v2": W / "templates_exhad"}


def link(dst: Path, src: Path):
    if not src.exists():
        raise FileNotFoundError(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.symlink_to(src.resolve())


def combined_v1(flavor, mass):
    parts = []
    for ch in CHANNELS:
        p = llp_csv_path(flavor, ch, mass, base=W / "llp_4vectors", mkdir=False)
        if not p.exists() or p.stat().st_size == 0:
            continue
        data = read_csv_matrix(p)[:, :5]
        if ch == "Bc":
            data = data.copy()
            data[:, 0] *= K_V1
        parts.append(data)
    return np.concatenate(parts)


for tag in TAGS:
    for flavor in FLAVORS:
        for mass in NEW + CHECK:
            label = format_mass_for_filename(mass)
            vec = W / f"vectors_{tag}" / flavor / "combined" / f"mN_{label}.csv"
            tpl = W / f"templates_{tag}" / flavor / f"templates_{label}.npz"
            if mass in CHECK:
                link(vec, llp_csv_path(flavor, "combined", mass, base=PUBLISHED[tag]["vectors"], mkdir=False))
                link(tpl, PUBLISHED[tag]["templates"] / flavor / f"templates_{label}.npz")
                continue
            produced = llp_csv_path(flavor, "combined", mass, base=W / "llp_4vectors", mkdir=False)
            if tag == "v2":
                link(vec, produced)
            else:
                vec.parent.mkdir(parents=True, exist_ok=True)
                write_csv_matrix(vec, combined_v1(flavor, mass))
            link(tpl, NEW_TEMPLATES[tag] / flavor / f"templates_{label}.npz")
print("assembled", len(NEW), "new and", len(CHECK), "check masses for", " ".join(TAGS))
