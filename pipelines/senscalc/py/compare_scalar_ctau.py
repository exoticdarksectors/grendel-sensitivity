"""Tracked BC4/BC5 scalar lifetime (grendel/models/scalar/data/eventcalc_1809/ctau.json.gz, exHad scalar-1809 =
EventCalc/SensCalc tables) against SensCalc's own ctauS.mx, decoded on this machine.

    python3 -I compare_scalar_ctau.py <ctauS_1.csv> <tracked ctau.json.gz>

SensCalc's evaluation (codes/Acceptances/Scalar.nb): log-log linear interpolation of the table. Reports the ratio
tracked/SensCalc at the tracked masses inside SensCalc's mass range: max |ratio - 1| overall and per mass decade.
"""
import bisect
import csv
import gzip
import json
import math
import sys


def main() -> int:
    rows = [r for r in csv.reader(open(sys.argv[1]))][1:]
    sm = [float(r[0]) for r in rows]
    sc = [float(r[1]) for r in rows]
    lm, lc = [math.log10(x) for x in sm], [math.log10(x) for x in sc]
    tracked = json.loads(gzip.decompress(open(sys.argv[2], "rb").read()))

    def senscalc(m):
        x = math.log10(m)
        i = min(max(bisect.bisect_right(lm, x) - 1, 0), len(lm) - 2)
        t = (x - lm[i]) / (lm[i + 1] - lm[i])
        return 10 ** (lc[i] + t * (lc[i + 1] - lc[i]))

    node = {round(m, 12): c for m, c in zip(sm, sc)}
    worst, n, nodes_hit, nodes_same = 0.0, 0, 0, 0
    per = {}
    for m, ct in tracked:
        if not (sm[0] <= m <= sm[-1]):
            continue
        r = ct / senscalc(m)
        n += 1
        d = abs(r - 1)
        worst = max(worst, d)
        dec = math.floor(math.log10(m))
        per[dec] = max(per.get(dec, 0.0), d)
        if round(m, 12) in node:
            nodes_hit += 1
            nodes_same += node[round(m, 12)] == ct
    print(f"SensCalc ctauS.mx: {len(sm)} nodes {sm[0]:g}-{sm[-1]:g} GeV; tracked: {len(tracked)} rows "
          f"{tracked[0][0]:g}-{tracked[-1][0]:g} GeV; compared {n} tracked masses inside the SensCalc range")
    print(f"tracked masses that are SensCalc nodes: {nodes_hit}, with identical ctau: {nodes_same}")
    print(f"max |ctau_tracked / ctau_SensCalc - 1| = {worst:.3e}")
    for dec in sorted(per):
        print(f"  m in [1e{dec}, 1e{dec + 1}) GeV: max {per[dec]:.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
