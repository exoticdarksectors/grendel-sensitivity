"""Progress of the running central grids (step 03) and a running comparison of the raw rows computed so far.

    source env/activate.sh && source pipelines/fonll/env.sh && python pipelines/fonll/progress.py

For each quark: raw quark-grid points done out of 100 x 80 (summed over the rapidity chunks), a linear ETA from
the chunk start times, and, for the points already written, the agreement with the stored macOS raw grid
(fonll-local/run/nlo_<q>/nl_<q>_.out, optional reference tree): byte-identical lines and the matched FONLL cross
section hdmv + G (hdrs - hdml) relative difference.
"""
import os
import time
from pathlib import Path

import numpy as np

RUN = Path(os.environ["FONLL_RUN_DIR"])
OUT = Path(os.environ["GRENDEL_FONLL_GRID_DIR"])
ORIG = Path(os.environ["FONLL_ORIG_WS"]) / "run"
NY, NPT = 100, 80
MASS = {"bottom": 4.75, "charm": 1.5}

for q in ("bottom", "charm"):
    final = OUT / f"fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_{q}.dat"
    if final.exists():
        print(f"{q}: finished ({final})")
        continue
    d = RUN / f"nlo_{q}_central"
    merged = d / f"nl_{q[0]}_.out"
    chunks = sorted(d.glob(f"grid_chunks*/*/nl_{q[0]}_.out"))
    if merged.exists() and len(merged.read_text().splitlines()) == 1 + NY * NPT:
        print(f"{q}: raw grid complete, fragmentation running")
        files = [merged]
    elif chunks:
        files = chunks
    else:
        print(f"{q}: not started")
        continue
    ref = ORIG / f"nlo_{q}" / f"nl_{q[0]}_.out"
    ref_rows = ref.read_text().splitlines()[1:] if ref.exists() else []  # the reference tree is optional
    ref_by_key = {tuple(l.split()[:2]): l for l in ref_rows}
    done, same, rel = 0, 0, []
    t_start = min(f.parent.stat().st_ctime for f in files)
    for f in files:
        for line in f.read_text().splitlines()[1:]:
            c = line.split()
            if len(c) < 10:
                continue
            done += 1
            r = ref_by_key.get(tuple(c[:2]))
            if r is None:
                continue
            same += line == r
            a = np.array([float(x) for x in c])
            b = np.array([float(x) for x in r.split()])
            g = b[0] ** 2 / (b[0] ** 2 + 25.0 * MASS[q] ** 2)
            xa, xb = a[2] + g * (a[4] - a[3]), b[2] + g * (b[4] - b[3])
            if xb != 0:
                rel.append(abs(xa - xb) / abs(xb))
    el = time.time() - t_start
    eta = el / done * (NY * NPT - done) if done else float("nan")
    rel = np.array(rel) if rel else np.array([0.0])
    cmp = (f"vs macOS: {same}/{done} lines identical, matched xsec |rel diff| max {rel.max():.2e}, "
           f"median {np.median(rel):.1e}" if ref_rows else "no stored macOS raw grid to compare with (GRENDEL_REF)")
    print(f"{q}: {done}/{NY * NPT} raw points ({100 * done / (NY * NPT):.1f} %), {el / 3600:.2f} h elapsed, "
          f"raw-grid ETA {eta / 3600:.2f} h (+ fragmentation); {cmp}")
