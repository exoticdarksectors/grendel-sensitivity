"""Exact point-in-fiducial test against the ray-cast mesh."""
from __future__ import annotations

from fractions import Fraction

import numpy as np

CHUNK = 20000
_ORIENT_ERRBOUND = (3.0 + 16.0 * 2.0 ** -53) * 2.0 ** -53


def _orient(ax, az, bx, bz, px, pz):
    """Exact sign (+1, 0, -1) of cross(b - a, p - a) in (x, z) for these doubles."""
    ax, az, bx, bz, px, pz = np.broadcast_arrays(ax, az, bx, bz, px, pz)
    t1 = (bx - ax) * (pz - az)
    t2 = (bz - az) * (px - ax)
    det = t1 - t2
    sign = np.sign(det)
    unsure = np.abs(det) <= _ORIENT_ERRBOUND * (np.abs(t1) + np.abs(t2))
    for k in (zip(*np.nonzero(unsure)) if unsure.any() else ()):
        fax, faz, fbx, fbz, fpx, fpz = (Fraction(float(c[k])) for c in (ax, az, bx, bz, px, pz))
        exact = (fbx - fax) * (fpz - faz) - (fbz - faz) * (fpx - fax)
        sign[k] = (exact > 0) - (exact < 0)
    return sign


def _side(ax, az, bx, bz, px, pz):
    """+1/-1: the side of the line a -> b the point lies on (``_orient``), a point exactly on the line
    taking the side of p + (eps, eps^2)."""
    s = _orient(ax, az, bx, bz, px, pz)
    zero = s == 0
    if zero.any():
        ax, az, bx, bz = (np.broadcast_to(c, s.shape)[zero] for c in (ax, az, bx, bz))
        s[zero] = np.where(bz != az, -np.sign(bz - az), np.sign(bx - ax))
    return s


class VerticalRayParity:
    """Point-in-mesh by the parity of the triangles above each point."""

    def __init__(self, mesh, cell=0.25):
        tri = np.asarray(mesh.triangles, dtype=np.float64)
        x, z = tri[:, :, 0], tri[:, :, 2]
        orient = _orient(x[:, 0], z[:, 0], x[:, 1], z[:, 1], x[:, 2], z[:, 2])
        keep = orient != 0
        tri = tri[keep]
        self.x0, self.y0, self.z0 = tri[:, :, 0], tri[:, :, 1], tri[:, :, 2]
        self.orient = orient[keep]
        self.area2 = ((self.x0[:, 1] - self.x0[:, 0]) * (self.z0[:, 2] - self.z0[:, 0])
                      - (self.x0[:, 2] - self.x0[:, 0]) * (self.z0[:, 1] - self.z0[:, 0]))
        self.cell = cell
        self.xmin = x.min() - 1.0
        self.zmin = z.min() - 1.0
        self.nx = int(np.ceil((x.max() + 1.0 - self.xmin) / cell)) + 1
        self.nz = int(np.ceil((z.max() + 1.0 - self.zmin) / cell)) + 1
        lists = [[] for _ in range(self.nx * self.nz)]
        bx0 = np.floor((self.x0.min(1) - self.xmin) / cell).astype(int)
        bx1 = np.floor((self.x0.max(1) - self.xmin) / cell).astype(int)
        bz0 = np.floor((self.z0.min(1) - self.zmin) / cell).astype(int)
        bz1 = np.floor((self.z0.max(1) - self.zmin) / cell).astype(int)
        for t in range(len(tri)):
            for ix in range(bx0[t], bx1[t] + 1):
                for iz in range(bz0[t], bz1[t] + 1):
                    lists[ix * self.nz + iz].append(t)
        kmax = max(len(cands) for cands in lists)
        self.cand = np.full((self.nx * self.nz, kmax), -1, dtype=np.int64)
        for c, cands in enumerate(lists):
            self.cand[c, :len(cands)] = cands
        self.ncand = np.array([len(cands) for cands in lists], dtype=np.int64)
        self.ymin, self.ymax = self.y0.min(axis=1), self.y0.max(axis=1)

    def crossings_above(self, points):
        pts = np.asarray(points, dtype=np.float64).reshape(-1, 3)
        out = np.zeros(len(pts), dtype=np.int64)
        ix = np.floor((pts[:, 0] - self.xmin) / self.cell).astype(np.int64)
        iz = np.floor((pts[:, 2] - self.zmin) / self.cell).astype(np.int64)
        idx = np.nonzero((ix >= 0) & (ix < self.nx) & (iz >= 0) & (iz < self.nz))[0]
        cell = ix[idx] * self.nz + iz[idx]
        order = np.argsort(self.ncand[cell], kind="stable")
        idx, cell = idx[order], cell[order]
        for s in range(0, len(idx), CHUNK):
            width = int(self.ncand[cell[s:s + CHUNK]].max())
            if width:
                out[idx[s:s + CHUNK]] = self._count(pts[idx[s:s + CHUNK]], self.cand[cell[s:s + CHUNK], :width])
        return out

    def _count(self, p, cand):
        """Crossings above each point ``p`` among its candidate triangles."""
        valid = cand >= 0
        ci = np.where(valid, cand, 0)
        X, Z = self.x0[ci], self.z0[ci]
        px, pz = p[:, 0][:, None], p[:, 2][:, None]
        orient = self.orient[ci]
        under = valid
        for i, j in ((0, 1), (1, 2), (2, 0)):
            under = under & (_side(X[..., i], Z[..., i], X[..., j], Z[..., j], px, pz) == orient)
        row, col = np.nonzero(under)
        t = ci[row, col]
        x, z, y = self.x0[t], self.z0[t], self.y0[t]
        qx, qz = p[row, 0], p[row, 2]
        with np.errstate(divide="ignore", invalid="ignore"):
            w1 = ((qx - x[:, 0]) * (z[:, 2] - z[:, 0]) - (x[:, 2] - x[:, 0]) * (qz - z[:, 0])) / self.area2[t]
            w2 = ((x[:, 1] - x[:, 0]) * (qz - z[:, 0]) - (qx - x[:, 0]) * (z[:, 1] - z[:, 0])) / self.area2[t]
            y_cross = (1.0 - w1 - w2) * y[:, 0] + w1 * y[:, 1] + w2 * y[:, 2]
        lo, hi = self.ymin[t], self.ymax[t]
        y_cross = np.where(np.isfinite(y_cross), np.clip(y_cross, lo, hi), 0.5 * (lo + hi))
        return np.bincount(row[y_cross > p[row, 1]], minlength=len(p))

    def __call__(self, points):
        return (self.crossings_above(points) % 2) == 1


_TEST = None


def points_in_fiducial(points):
    """Bool array: True where each (M, 3) point lies inside the fiducial mesh (the ray-cast one),
    exactly; drop-in for the analytic test in ``grendel_geometry``."""
    global _TEST
    if _TEST is None:
        from .grendel_geometry import mesh_fiducial
        _TEST = VerticalRayParity(mesh_fiducial)
    return _TEST(points)
