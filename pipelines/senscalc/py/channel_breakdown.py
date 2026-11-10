"""Per-channel yields of a SensCalc 'Tabulated Nevents' folder at chosen masses (log-log interpolation as in SensCalc).
python3 -I channel_breakdown.py <Nevents folder> <threshold> m1 m2 ...  (baseline channels only: AP-Central, DY-Central)"""
import sys, glob, math, collections
folder, thr = sys.argv[1], float(sys.argv[2]); masses = [float(x) for x in sys.argv[3:]]
tabs = {}
for f in sorted(glob.glob(folder + '/Nevents_*.dat')):
    prod = f.split('/')[-1].split('_')[2]
    if prod.endswith('Lower') or prod.endswith('Upper'):
        continue
    d = collections.defaultdict(dict)
    for line in open(f):
        p = line.split()
        if len(p) == 3:
            d[float(p[0])][float(p[1])] = float(p[2])
    tabs[prod] = d
prods = sorted(tabs)
def nev(prod, m, y):
    d = tabs[prod]; ms = sorted(d)
    if not (ms[0] <= m <= ms[-1]): return 0.0
    i = max(k for k in range(len(ms)) if ms[k] <= m); j = min(i + 1, len(ms) - 1)
    def at(mm):
        ys = sorted(d[mm])
        if not (ys[0] <= y <= ys[-1]): return 0.0
        a = max(k for k in range(len(ys)) if ys[k] <= y); b = min(a + 1, len(ys) - 1)
        la, lb = math.log10(d[mm][ys[a]] + 1e-90), math.log10(d[mm][ys[b]] + 1e-90)
        t = 0 if a == b else (math.log10(y) - math.log10(ys[a])) / (math.log10(ys[b]) - math.log10(ys[a]))
        return 10 ** (la + t * (lb - la))
    if i == j: return at(ms[i])
    t = (math.log10(m) - math.log10(ms[i])) / (math.log10(ms[j]) - math.log10(ms[i]))
    return 10 ** ((1 - t) * math.log10(at(ms[i]) + 1e-90) + t * math.log10(at(ms[j]) + 1e-90))
ygrid = [10 ** (-11 + 0.017 * k) for k in range(589)]
print('channels:', prods)
for m in masses:
    tot = [sum(nev(p, m, y) for p in prods) for y in ygrid]
    above = [k for k, t in enumerate(tot) if t >= thr]
    if not above:
        print('m=%.2f: no sensitivity, max N_tot = %.3g' % (m, max(tot))); continue
    k = above[0]; y = ygrid[k]
    shares = {p: nev(p, m, y) / tot[k] for p in prods}
    peaks = {p: max(nev(p, m, yy) for yy in ygrid[::3]) for p in prods}
    print('m=%.2f: lower edge g_Y~%.3g, N_tot max %.3g; shares at edge: %s' % (m, y, max(tot),
          ', '.join('%s %.2f' % (p, s) for p, s in sorted(shares.items(), key=lambda x: -x[1]) if s > 0.01)))
    print('        peak N per channel: ' + ', '.join('%s %.3g' % (p, v) for p, v in sorted(peaks.items(), key=lambda x: -x[1]) if v > 0.01))
