"""Smoke-test criterion for the headless SensCalc chain (10_smoke.sh): does the chain run and land in the right place?

    python3 -I smoke_check.py <smoke tree> <campaign tree> <experiment>

The smoke run draws the notebook-2 meson-decay spectra with 2e4 instead of 1e6 events per channel. FASER2-FPF sees only
the very forward tail of the B -> K(*) a sample (theta < 1.5 mrad: 0.15 % of it; its long-lifetime edge mostly the
E < 300 GeV part, 0.02-0.04 %), i.e. a few to a few tens of sampled ALPs per mass at 2e4: the B-channel N_ev at the
FASER2 edges then fluctuates by factors of several (observed 0.004x-15x), so the curve itself cannot be held to a
percent-level tolerance at this statistic (a 5 % criterion failed for that reason alone). The criterion therefore is:

1. Deterministic channels, exact (catches environment bugs such as the case-sensitive path of 2026-11-07, which made
   the bremsstrahlung N_ev 0.001-0.03x the campaign's):
     bremsstrahlung (AP-Lower/Central/Upper): N_ev ratio to the campaign within 0.5 % at every (m, y) with N_ref >= 1;
     fragmentation (pregenerated spectrum):   within 2 %.
2. Placement of the curve, sized for the statistic: both edges of the central variant within 50 % of the campaign at the
   17 masses of compare_curves where both have sensitivity, and the mass reach within 5 %.
Exit status 1 if any criterion fails; every number is printed.
"""
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import contour_tools as ct  # noqa: E402
from compare_curves import MASSES  # noqa: E402


def load_nev(f):
    d = {}
    for line in open(f):
        p = line.split()
        if len(p) == 3:
            try:
                d[(round(float(p[0]), 8), round(float(p[1]), 12))] = float(p[2])
            except ValueError:
                pass
    return d


def main() -> int:
    smoke, camp, exp = sys.argv[1], sys.argv[2], sys.argv[3]
    ok = True
    nev = "Tabulated Nevents/ALP-fermion/" + exp
    for pattern, tol, label in (("Bremsstrahlung-AP-*", 0.005, "bremsstrahlung"), ("Fragmentation", 0.02, "fragmentation")):
        for f in sorted(glob.glob(os.path.join(smoke, nev, f"Nevents_ALP-fermion_{pattern}_at_*.dat"))):
            g = os.path.join(camp, nev, os.path.basename(f))
            a, b = load_nev(f), load_nev(g)
            r = [a[k] / b[k] for k in b if k in a and b[k] >= 1.0]
            worst = max((abs(x - 1) for x in r), default=float("nan"))
            good = bool(r) and worst <= tol
            ok &= good
            print(f"{label:15s} {os.path.basename(f).split('_')[2]:28s} points {len(r):6d}  max |ratio-1| {100 * worst:.3f}%"
                  f"  (tol {100 * tol:.1f}%)  {'PASS' if good else 'FAIL'}")
    sd = "Sensitivity domains/ALP-fermion/" + exp
    for f in sorted(glob.glob(os.path.join(smoke, sd, "Sensitivity_*.json"))):
        if f.endswith(".metadata.json"):
            continue
        t = ct.load_senscalc_json(f)[0]
        r = ct.load_senscalc_json(os.path.join(camp, sd, os.path.basename(f)))[0]
        devs = []
        for m in MASSES:
            et, er = ct.edges(t, m), ct.edges(r, m)
            if et and er:
                devs += [(et[0] / er[0] - 1, m, "lower"), (et[1] / er[1] - 1, m, "upper")]
        worst = max(devs, key=lambda x: abs(x[0]))
        reach_t, reach_r = ct.mass_range(t)[1], ct.mass_range(r)[1]
        good = abs(worst[0]) <= 0.5 and abs(reach_t / reach_r - 1) <= 0.05
        ok &= good
        print(f"curve {os.path.basename(f)}: max edge deviation {100 * worst[0]:+.1f}% ({worst[2]}, m={worst[1]}) (tol 50%), "
              f"reach {reach_t:.3f} vs {reach_r:.3f} GeV ({100 * (reach_t / reach_r - 1):+.1f}%, tol 5%)  {'PASS' if good else 'FAIL'}")
    print("SMOKE:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
