"""Helpers for BC10 sensitivity contours (closed chains in the (mass, coupling) plane).

Coupling conventions: SensCalc ALP-fermion y = g_Y = 2 v_h / f_a (arXiv:2501.04525); the GRENDEL curve files store
1/f in the BNT convention = g_Y / 246 (header of ship_senscalc.dat).
"""
import json
import math

GY_TO_INVF = 1.0 / 246.0


def load_dat(path):
    """Chains from a 'mass  coupling' file; chains separated by blank lines, '#' comments ignored."""
    chains, cur = [], []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith('#'):
                continue
            if not s:
                if cur:
                    chains.append(cur)
                    cur = []
                continue
            a, b = s.split()[:2]
            cur.append((float(a), float(b)))
    if cur:
        chains.append(cur)
    return chains


def load_senscalc_json(path, brem='AP-Central', dy='Central', decay='Default'):
    """Chains (mass, g_Y) of one entry of a SensCalc Sensitivity_*.json file."""
    with open(path) as fh:
        data = json.load(fh)
    for entry in data:
        if (entry.get('Bremsstrahlung description') == brem and entry.get('Drell-Yan description') == dy
                and entry.get('Decay description') == decay):
            return [[(p[0], p[1]) for p in ch] for ch in entry['Sensitivity domains']], entry['Production modes']
    raise KeyError((brem, dy, decay))


def crossings(chains, m):
    """Coupling values where the closed chains cross the vertical line at mass m (log-log interpolation)."""
    ys = []
    lm = math.log(m)
    for ch in chains:
        pts = list(ch)
        if pts[0] != pts[-1]:
            pts.append(pts[0])
        for (x1, y1), (x2, y2) in zip(pts[:-1], pts[1:]):
            if x1 <= 0 or x2 <= 0 or y1 <= 0 or y2 <= 0:
                continue
            l1, l2 = math.log(x1), math.log(x2)
            if l1 == l2:
                continue
            if (l1 - lm) * (l2 - lm) <= 0 and not (l2 == lm and l1 != lm and False):
                t = (lm - l1) / (l2 - l1)
                if 0 <= t <= 1:
                    ys.append(math.exp(math.log(y1) + t * (math.log(y2) - math.log(y1))))
    return sorted(set(round(y, 300) for y in ys))


def edges(chains, m):
    """(lower, upper) edge at mass m, or None if the line does not cross the region."""
    ys = crossings(chains, m)
    if len(ys) < 2:
        return None
    return ys[0], ys[-1]


def mass_range(chains):
    xs = [p[0] for ch in chains for p in ch]
    return min(xs), max(xs)


def scale(chains, f):
    return [[(x, y * f) for x, y in ch] for ch in chains]


def write_dat(path, chains, header_lines):
    with open(path, 'w') as fh:
        for h in header_lines:
            fh.write('# ' + h + '\n')
        for i, ch in enumerate(chains):
            if i:
                fh.write('\n')
            fh.write('# chain: %d vertices\n' % len(ch))
            for x, y in ch:
                fh.write('%.6e  %.6e\n' % (x, y))
