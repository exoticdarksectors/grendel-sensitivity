"""Vector extraction from PDF path operators (PyMuPDF).

A figure's curves are stroke (``type 's'``) and fill (``'f'``/``'fs'``)
drawings; each is a list of items ``('l', p1, p2)``, ``('m', p)``,
``('c', p0, p1, p2, p3)`` or ``('re', rect)``. The helpers here flatten those
into point chains, select drawings by colour, join chains whose endpoints
touch, and write the chains out in the blank-line-separated ``.dat``
format the renderer reads.

Bezier segments are sampled at the parameter values in ``t_samples``. The
extractors pin their own sampling (four quarter points, or ``linspace`` of
6, 7 or 8) because the tracked files were produced with it.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

QUARTERS = (0.25, 0.5, 0.75, 1.0)


def linspace_samples(n: int) -> tuple[float, ...]:
    return tuple(float(t) for t in np.linspace(1.0 / n, 1.0, n))


def bezier(p0, p1, p2, p3, t):
    x = (1 - t) ** 3 * p0.x + 3 * (1 - t) ** 2 * t * p1.x + 3 * (1 - t) * t * t * p2.x + t ** 3 * p3.x
    y = (1 - t) ** 3 * p0.y + 3 * (1 - t) ** 2 * t * p1.y + 3 * (1 - t) * t * t * p2.y + t ** 3 * p3.y
    return x, y


def near(a, b, tol=1e-3):
    return all(abs(x - y) < tol for x, y in zip(a, b))


def color_matches(color, target, tol=5e-3):
    col = tuple(round(c, 3) for c in color)
    return max(abs(a - b) for a, b in zip(col, target)) <= tol


def is_dashed(drawing) -> bool:
    return drawing.get("dashes") not in (None, "", "[] 0")


def flatten(items, t_samples=QUARTERS, *, rects=False):
    """One chain per drawing: a segment whose start is not the running pen
    position starts from its own first point."""
    pts = []
    for it in items:
        tag = it[0]
        if tag == "l":
            p1, p2 = it[1], it[2]
            if not pts or pts[-1] != (p1.x, p1.y):
                pts.append((p1.x, p1.y))
            pts.append((p2.x, p2.y))
        elif tag == "m":
            pts.append((it[1].x, it[1].y))
        elif tag == "c":
            p0 = it[1]
            if not pts or pts[-1] != (p0.x, p0.y):
                pts.append((p0.x, p0.y))
            for t in t_samples:
                pts.append(bezier(it[1], it[2], it[3], it[4], t))
        elif tag == "re" and rects:
            r = it[1]
            pts.extend([(r.x0, r.y0), (r.x1, r.y0), (r.x1, r.y1), (r.x0, r.y1), (r.x0, r.y0)])
    return pts


def stitch(items, t_samples=QUARTERS):
    """Like ``flatten`` for a stroke known to be one continuous polyline:
    line items contribute their end point only."""
    pts = []
    for it in items:
        tag = it[0]
        if tag == "l":
            p1, p2 = it[1], it[2]
            if not pts:
                pts.append((p1.x, p1.y))
            pts.append((p2.x, p2.y))
        elif tag == "m":
            pts.append((it[1].x, it[1].y))
        elif tag == "c":
            p0 = it[1]
            if not pts:
                pts.append((p0.x, p0.y))
            for t in t_samples:
                pts.append(bezier(it[1], it[2], it[3], it[4], t))
    return pts


def subpaths(items, t_samples=QUARTERS):
    """Split a drawing into subpaths at every ``m`` and at every segment
    that does not continue the running pen position (a union fill of
    several polygons)."""
    out = []
    cur = []
    last = None

    def flush():
        nonlocal cur
        if cur:
            out.append(cur)
            cur = []

    for it in items:
        tag = it[0]
        if tag == "l":
            p1, p2 = it[1], it[2]
            if not cur or last != (p1.x, p1.y):
                flush()
                cur.append((p1.x, p1.y))
            cur.append((p2.x, p2.y))
            last = (p2.x, p2.y)
        elif tag == "m":
            flush()
            cur.append((it[1].x, it[1].y))
            last = (it[1].x, it[1].y)
        elif tag == "c":
            p0, p3 = it[1], it[4]
            if not cur or last != (p0.x, p0.y):
                flush()
                cur.append((p0.x, p0.y))
            for t in t_samples:
                cur.append(bezier(it[1], it[2], it[3], it[4], t))
            last = (p3.x, p3.y)
    flush()
    return out


def collect_strokes(page, color, *, t_samples, tol=5e-3, dashed=None, width=None, opaque_only=False,
                    min_points=2):
    """Every stroked chain of one colour. ``dashed`` selects solid (False)
    or dashed (True) strokes; ``opaque_only`` drops translucent band
    edges; ``width`` pins the line width."""
    chains = []
    for d in page.get_drawings():
        if d["type"] != "s" or not d["color"]:
            continue
        if not color_matches(d["color"], color, tol):
            continue
        if dashed is not None and is_dashed(d) != dashed:
            continue
        if width is not None and abs((d.get("width") or 0) - width) > 0.2:
            continue
        if opaque_only and (d.get("stroke_opacity") or 1.0) < 0.99:
            continue
        pts = flatten(d["items"], t_samples)
        if len(pts) >= min_points:
            chains.append(pts)
    return chains


def collect_fills(page, color, *, t_samples, opacity=None, tol=5e-3, min_points=3, rects=False):
    polys = []
    for d in page.get_drawings():
        if not d["fill"]:
            continue
        if not color_matches(d["fill"], color, tol):
            continue
        if opacity is not None and abs((d.get("fill_opacity") or 1.0) - opacity) > 0.05:
            continue
        pts = flatten(d["items"], t_samples, rects=rects)
        if len(pts) >= min_points:
            polys.append(pts)
    return polys


def merge_chains(chains, tol=0.75):
    """Join polylines whose endpoints (nearly) touch."""
    chains = [list(c) for c in chains]

    def d2(p, q):
        return (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2

    merged = True
    while merged:
        merged = False
        for i in range(len(chains)):
            for j in range(len(chains)):
                if i == j:
                    continue
                a, b = chains[i], chains[j]
                if d2(a[-1], b[0]) < tol ** 2:
                    a.extend(b[1:])
                elif d2(a[-1], b[-1]) < tol ** 2:
                    a.extend(reversed(b[:-1]))
                elif d2(a[0], b[-1]) < tol ** 2:
                    chains[i] = b + a[1:]
                elif d2(a[0], b[0]) < tol ** 2:
                    chains[i] = list(reversed(b)) + a[1:]
                else:
                    continue
                chains.pop(j)
                merged = True
                break
            if merged:
                break
    return chains


def write_chains(path, chains, header_lines, *, precision=6, label="chain", numbered=False, columns=None):
    """Blank-line-separated chains of already converted ``(x, y)`` rows.

    ``columns`` is the comment naming the two columns; ``numbered`` writes
    ``# chain k: N vertices`` instead of ``# chain: N vertices``.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fmt = f"{{:.{precision}e}}  {{:.{precision}e}}\n"
    with path.open("w") as fh:
        for line in header_lines:
            fh.write(f"# {line}\n")
        if columns:
            fh.write(f"# {columns}\n")
        for k, chain in enumerate(chains):
            if k:
                fh.write("\n")
            fh.write(f"# {label} {k}: {len(chain)} vertices\n" if numbered else f"# {label}: {len(chain)} vertices\n")
            for x, y in chain:
                fh.write(fmt.format(x, y))
    print(f"wrote {path} ({len(chains)} chains)")


def clip_chains(chains, to_data, clip_x=None):
    """Convert each chain to data space, splitting it into contiguous
    in-frame runs at clip crossings rather than dropping vertices in place
    (deleting an off-frame arc would let the plotter bridge its neighbours
    with a spurious chord)."""
    out = []
    for ch in chains:
        run = []
        for x, y in ch:
            m, v = to_data(x, y)
            if clip_x and not (clip_x[0] <= m <= clip_x[1]):
                if run:
                    out.append(run)
                    run = []
                continue
            run.append((m, v))
        if run:
            out.append(run)
    return out
