"""Input trees for the 4x-statistics scans of 2.6-3.7 GeV: the vectors produced here
(v2 as produced, v1 with the B_c rows at the v1 rate 0.9e6/1.64e6) and the existing
templates (published set where it has the mass, else the 2026-10-06 refinement's)."""
import sys
from pathlib import Path
import numpy as np
from grendel.io.vectors import format_mass_for_filename
from grendel.models.hnl.production.combine_channels import CHANNELS
from grendel.models.hnl.production.paths import llp_csv_path, read_csv_matrix, write_csv_matrix

H = Path(__file__).resolve().parent
C = H.parents[1]
MW = [float(m) for m in sys.argv[1].split()]
TEMPLATES = {"v1": [C / "hnl/hnl_exhad_20260917/tmp_fairship300/decay_templates", C / "hnl/hnl_refine_20261006/templates_fairship"],
             "v2": [C / "v2_baseline_20260929/templates/hnl", C / "hnl/hnl_refine_20261006/templates_exhad"]}

def link(dst, src):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.symlink_to(src)

for flavor in ("Ue", "Umu", "Utau"):
    for m in MW:
        label = format_mass_for_filename(m)
        produced = llp_csv_path(flavor, "combined", m, base=H / "llp_4vectors", mkdir=False)
        link(H / "vectors_v2" / flavor / "combined" / produced.name, produced)
        parts = []
        for ch in CHANNELS:
            p = llp_csv_path(flavor, ch, m, base=H / "llp_4vectors", mkdir=False)
            if p.exists() and p.stat().st_size > 0:
                d = read_csv_matrix(p)[:, :5]
                if ch == "Bc":
                    d = d.copy(); d[:, 0] *= 0.9e6 / 1.64e6
                parts.append(d)
        v1 = H / "vectors_v1" / flavor / "combined" / produced.name
        v1.parent.mkdir(parents=True, exist_ok=True)
        write_csv_matrix(v1, np.concatenate(parts))
        for tag, dirs in TEMPLATES.items():
            src = next(d / flavor / f"templates_{label}.npz" for d in dirs if (d / flavor / f"templates_{label}.npz").exists())
            link(H / f"templates_{tag}" / flavor / f"templates_{label}.npz", src)
print("assembled", len(MW), "masses x 3 flavours for v1 and v2")
