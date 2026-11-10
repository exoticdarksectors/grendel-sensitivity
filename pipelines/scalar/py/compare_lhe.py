"""Compare two LHE files event by event (same MG5 seed): count, and the largest relative difference of the
four-momenta of every particle line.

    python compare_lhe.py OURS.lhe.gz REF.lhe.gz
"""
import gzip
import sys

import numpy as np


def events(path):
    with gzip.open(path, "rt") as fh:
        inside, ev = False, []
        for line in fh:
            if line.startswith("<event"):
                inside, ev = True, []
                continue
            if line.startswith("</event"):
                inside = False
                yield ev
                continue
            if inside:
                parts = line.split()
                if len(parts) == 13:            # particle line: id st m1 m2 c1 c2 px py pz E m vtim spin
                    ev.append([float(x) for x in parts[6:10]] + [int(parts[0])])


n = n_diff = n_ids = 0
worst = 0.0
for a, b in zip(events(sys.argv[1]), events(sys.argv[2])):
    n += 1
    A, B = np.array(a), np.array(b)
    if A.shape != B.shape or not np.array_equal(A[:, 4], B[:, 4]):
        n_ids += 1
        continue
    if not np.array_equal(A, B):
        n_diff += 1
        scale = np.maximum(np.abs(B[:, 3:4]), 1e-9)          # relative to each particle's energy
        worst = max(worst, float(np.max(np.abs(A[:, :4] - B[:, :4]) / scale)))
print(f"{n} events compared: {n - n_diff - n_ids} identical, {n_diff} with the same particles and momenta "
      f"differing by at most {worst:.2g} of the particle energy, {n_ids} with different particle content")
