"""Cross-check of the fragmentation half of the chain: the stored macOS raw quark grids, fragmented here.

    source env/activate.sh && source pipelines/fonll/env.sh
    python pipelines/common/slot_run.py -n 1 --label fonll:xcheck -- python pipelines/fonll/xcheck_fragment.py

Copies the raw quark grids the tracked grids were made from (fonll-local/run/nlo_{bottom,charm}/nl_{b,c}_.out,
2026-05-31, byte-identical to the 2026-06-10 reruns) and the original charm calibration
(fonll-local/output/charm_feeddown_calibration.json, weight 1.2152533641182024) into work/fonll/xcheck/, then
runs the shipped generate() with reuse_existing_grids=True, i.e. this machine's fragmfonll, the D* feeddown
combination and the grid writer. With identical raw grids the result isolates the platform dependence of
everything after fonllgridlha; compared with the tracked grids by 05_compare.py (section "xcheck") and here.
Minutes of CPU. Needs the optional reference tree (GRENDEL_REF); exits after a SKIP line without it.
"""
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np

from grendel.production.fonll_grids import generate as gen

sys.path.insert(0, str(Path(os.environ["REPRO"]) / "pipelines/common"))
from ref import skip_unless  # noqa: E402

FONLL_WORK = Path(os.environ["FONLL_WORK"])
ORIG = Path(os.environ["FONLL_ORIG_WS"])
TRACKED = Path(os.environ["FONLL_TRACKED_DIR"])
skip_unless("FONLL xcheck: the stored macOS raw grids fragmented here", ORIG / "output/charm_feeddown_calibration.json",
            ORIG / "run/nlo_bottom/nl_b_.out", ORIG / "run/nlo_charm/nl_c_.out")
root = FONLL_WORK / "xcheck"
ws = gen.Workspace.resolve(None, root / "output", root / "run", root / "logs")
ws.out.mkdir(parents=True, exist_ok=True)
shutil.copy2(ORIG / "output" / "charm_feeddown_calibration.json", ws.out / "charm_feeddown_calibration.json")
calib = json.loads((ws.out / "charm_feeddown_calibration.json").read_text())

for q in ("bottom", "charm"):
    run_dir = ws.run / f"nlo_{q}_central"
    run_dir.mkdir(parents=True, exist_ok=True)
    raw = run_dir / f"nl_{q[0]}_.out"
    shutil.copy2(ORIG / "run" / f"nlo_{q}" / raw.name, raw)
    s = gen.generate(ws, "nlo", q, gen.central_variation(q), reuse_existing_grids=True, grid_workers=1,
                     feeddown_calibration=calib if q == "charm" else None)
    assert s["raw_grid_reused"], "the stored raw grid was not reused"
    name = gen.grid_filename("nlo", "central", q)
    a = np.loadtxt(ws.out / name, comments="#")
    b = np.loadtxt(TRACKED / name, comments="#")
    same_rows = sum(x == y for x, y in zip([l for l in (ws.out / name).read_text().splitlines() if not l.startswith("#")],
                                           [l for l in (TRACKED / name).read_text().splitlines() if not l.startswith("#")]))
    nz = b[:, 2] != 0
    rel = np.abs(a[nz, 2] - b[nz, 2]) / np.abs(b[nz, 2])
    print(f"[xcheck {q}] stored macOS raw grid fragmented here vs tracked grid: {same_rows}/10000 data rows "
          f"byte-identical; max |rel diff| {rel.max():.2e}, median {np.median(rel):.1e}")
