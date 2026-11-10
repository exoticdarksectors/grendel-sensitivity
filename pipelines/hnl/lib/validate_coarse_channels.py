"""Compare the coarse-grid channel files of step 05/06 with the stored published inputs.

Dmeson, Bc, Kmeson (as-built): file and row by row against asbuilt_20261014/hnl/llp_4vectors/<fl>/<ch>/.
Bmeson, Bbaryon, induced_tau (June v1, deterministic): row by row against their blocks in the stored combined files
  (rows identical as text; rows whose kinematics differ, and by how much; weight sums).
tau (prompt; MadGraph pool) and WZ (MadGraph): weight sums and mean energies against the stored blocks.
Writes work/hnl/validation/coarse_channels.csv and prints a summary.   python validate_coarse_channels.py [masses...]
The stored four-vectors are in the optional reference tree ($GRENDEL_REF): without them this prints a SKIP line.
"""
import hashlib
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPRO = Path(os.environ["REPRO"])
sys.path.insert(0, str(REPRO / "pipelines/hnl/lib"))
import grids  # noqa: E402
from stored import ASBUILT_VECTORS, read_rows, split_stored_combined, to_array  # noqa: E402
from ref import skip_unless  # noqa: E402  (pipelines/common, on sys.path through stored)

skip_unless("HNL coarse channels vs the stored as-built and June four-vectors", ASBUILT_VECTORS)
C = REPRO / "work/hnl/coarse"
masses = [float(x) for x in sys.argv[1:]] or grids.COARSE
rows = []
for fl in grids.FLAVORS:
    for m in masses:
        lab = grids.label(m)
        mine = {}
        for prod, chans in (("asbuilt", ["Dmeson", "Bc", "Kmeson"]), ("v1", ["Bmeson", "Bbaryon", "tau", "induced_tau"]),
                            ("wz", ["WZ"])):
            for ch in chans:
                mine[ch] = C / prod / "llp_4vectors" / fl / ch / f"mN_{lab}.csv"
        if not all(p.exists() for p in mine.values()):
            continue
        for ch in ("Dmeson", "Bc", "Kmeson"):
            r_mine, r_ref = read_rows(mine[ch]), read_rows(ASBUILT_VECTORS / fl / ch / f"mN_{lab}.csv")
            A, B = to_array(r_mine), to_array(r_ref)
            same_n = len(A) == len(B)
            rows.append(dict(flavor=fl, mass=m, channel=ch, kind="rows", n=len(A), n_ref=len(B), identical=r_mine == r_ref,
                             wsum=A[:, 0].sum() if len(A) else 0.0, wsum_ref=B[:, 0].sum() if len(B) else 0.0,
                             rows_text_identical=sum(x == y for x, y in zip(r_mine, r_ref)) if same_n else -1,
                             rows_kin_differ=int(np.any(A[:, 1:] != B[:, 1:], axis=1).sum()) if same_n else -1,
                             max_rel_kin=float(np.max(np.abs(A - B) / np.maximum(np.abs(B), 1e-300)))
                             if same_n and len(A) else 0.0))
        blocks = split_stored_combined(fl, m, len(read_rows(mine["Bmeson"])), len(read_rows(mine["Bbaryon"])),
                                       n_wz=len(read_rows(mine["WZ"])) and 100_000)
        for ch in ("Bmeson", "Bbaryon", "induced_tau", "tau", "WZ"):
            r_mine, r_ref = read_rows(mine[ch]), blocks[ch]
            A, B = to_array(r_mine), to_array(r_ref)
            d = dict(flavor=fl, mass=m, channel=ch, n=len(A), n_ref=len(B),
                     wsum=A[:, 0].sum() if len(A) else 0.0, wsum_ref=B[:, 0].sum() if len(B) else 0.0)
            same_n = len(A) == len(B)
            if ch in ("tau", "WZ"):
                d.update(kind="stat", meanE=A[:, 1].mean() if len(A) else np.nan,
                         meanE_ref=B[:, 1].mean() if len(B) else np.nan, identical=r_mine == r_ref,
                         rows_text_identical=sum(x == y for x, y in zip(r_mine, r_ref)) if same_n else -1)
            else:
                d.update(kind="rows", identical=r_mine == r_ref,
                         rows_text_identical=sum(x == y for x, y in zip(r_mine, r_ref)) if same_n else -1,
                         rows_kin_differ=int(np.any(A[:, 1:] != B[:, 1:], axis=1).sum()) if same_n else -1,
                         max_rel_kin=float(np.max(np.abs(A[:, 1:] - B[:, 1:]) / np.maximum(np.abs(B[:, 1:]), 1e-300)))
                         if same_n and len(A) else 0.0)
            rows.append(d)
t = pd.DataFrame(rows)
out = REPRO / "work/hnl/validation"
out.mkdir(parents=True, exist_ok=True)
t.to_csv(out / "coarse_channels.csv", index=False)
for ch, g in t.groupby("channel", sort=False):
    if g.kind.iloc[0] == "bytes":
        print(f"{ch:12s} byte-identical {int(g.identical.sum())}/{len(g)} files")
    elif g.kind.iloc[0] == "rows":
        nz = g[g.n_ref > 0]
        print(f"{ch:12s} identical files {int(g.identical.sum())}/{len(g)}; rows text-identical "
              f"{int(nz.rows_text_identical.clip(lower=0).sum())}/{int(nz.n_ref.sum())}; rows with kinematics "
              f"differing {int(nz.rows_kin_differ.clip(lower=0).sum())} (max rel diff {nz.max_rel_kin.max():.1e}); "
              f"w_sum ratio {np.nanmin(nz.wsum / nz.wsum_ref):.7f}-{np.nanmax(nz.wsum / nz.wsum_ref):.7f}")
    else:
        nz = g[g.n_ref > 0]
        r = nz.wsum / nz.wsum_ref
        print(f"{ch:12s} {len(nz)} non-empty points: identical files {int(nz.identical.sum())}, rows text-identical "
              f"{int(nz.rows_text_identical.clip(lower=0).sum())}/{int(nz.n_ref.sum())}; w_sum ratio median "
              f"{r.median():.6f} ({r.min():.6f}-{r.max():.6f}); <E> ratio median {(nz.meanE / nz.meanE_ref).median():.5f}")
