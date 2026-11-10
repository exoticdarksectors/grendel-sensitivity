"""Compare SensCalc ALP-fermion sensitivity regions (Sensitivity_*.json) edge by edge.

    python3 -I compare_curves.py <test.json|.dat> <reference.json|.dat> [--label NAME] [--out file.tsv] [--masses m1,m2,...]
                                 [--tol 0.01] [--variants all|central]

For every variant (Bremsstrahlung x Drell-Yan description; "central" = AP-Central + Central) present in both files,
the lower and upper edge in g_Y at each mass (closed domains, log-log interpolation along the boundary, the
contour_tools.edges of the campaign) and the mass range are compared. Prints the largest |test/ref - 1| per variant;
exit status 1 if the central variant deviates by more than --tol at any mass where both have sensitivity, or if the
sets of masses with sensitivity differ.
"""
import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import contour_tools as ct  # noqa: E402  (campaign helper, copied verbatim)

# 17 masses, 0.05-4.5 GeV, away from SensCalc's resonance windows (0.125-0.14, 0.538-0.555, 0.94-0.974 GeV)
MASSES = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
KEYS = ("Bremsstrahlung description", "Drell-Yan description", "Decay description")


def variants(path):
    if path.endswith(".dat"):  # GRENDEL curve file (chains of mass, 1/f): one variant, the central one
        return {("AP-Central", "Central", "Default"): ct.load_dat(path)}
    out = {}
    for e in json.load(open(path)):
        key = tuple(e.get(k) for k in KEYS)
        out[key] = [[(p[0], p[1]) for p in ch] for ch in e["Sensitivity domains"]]
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("test")
    ap.add_argument("ref")
    ap.add_argument("--label", default="")
    ap.add_argument("--out")
    ap.add_argument("--masses")
    ap.add_argument("--tol", type=float, default=0.01)
    ap.add_argument("--variants", default="all", choices=["all", "central"])
    a = ap.parse_args()
    masses = [float(x) for x in a.masses.split(",")] if a.masses else MASSES
    t, r = variants(a.test), variants(a.ref)
    common = [k for k in r if k in t]
    if a.variants == "central":
        common = [k for k in common if k[0] == "AP-Central" and k[1] == "Central"]
    rows = ["variant\tm_GeV\tlower_ref\tlower_test/ref\tupper_ref\tupper_test/ref"]
    status = 0
    print(f"{a.label}: test {os.path.basename(a.test)}  ref {os.path.basename(a.ref)}  variants: test {len(t)}, ref {len(r)}, "
          f"common {len(common)}")
    for k in common:
        name = f"{k[0]}|{k[1]}"
        worst_lo = worst_up = 0.0
        where_lo = where_up = None
        mismatch = []
        for m in masses:
            et, er = ct.edges(t[k], m), ct.edges(r[k], m)
            if (et is None) != (er is None):
                mismatch.append(m)
                rows.append(f"{name}\t{m}\t{er[0] if er else '-'}\t-\t{er[1] if er else '-'}\t-")
                continue
            if et is None:
                rows.append(f"{name}\t{m}\t-\t-\t-\t-")
                continue
            rl, ru = et[0] / er[0], et[1] / er[1]
            rows.append(f"{name}\t{m}\t{er[0]:.5e}\t{rl:.5f}\t{er[1]:.5e}\t{ru:.5f}")
            if abs(rl - 1) > abs(worst_lo):
                worst_lo, where_lo = rl - 1, m
            if abs(ru - 1) > abs(worst_up):
                worst_up, where_up = ru - 1, m
        (lt, ht), (lr, hr) = ct.mass_range(t[k]), ct.mass_range(r[k])
        rows.append(f"{name}\trange\t{lr:.4f}-{hr:.4f}\t{lt:.4f}-{ht:.4f}\t\t")
        central = k[0] == "AP-Central" and k[1] == "Central"
        print(f"  {name:28s} max|lower-1| {100 * worst_lo:+.2f}% (m={where_lo})  max|upper-1| {100 * worst_up:+.2f}% "
              f"(m={where_up})  range test {lt:.4f}-{ht:.4f} ref {lr:.4f}-{hr:.4f}"
              + (f"  SENSITIVITY MISMATCH at {mismatch}" if mismatch else "") + ("  <- central" if central else ""))
        if central and (mismatch or max(abs(worst_lo), abs(worst_up)) > a.tol):
            status = 1
    if a.out:
        os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
        open(a.out, "w").write("\n".join(rows) + "\n")
    print(f"  central within {100 * a.tol:.2f}%: {'PASS' if status == 0 else 'FAIL'}")
    return status


if __name__ == "__main__":
    sys.exit(main())
