"""Figs. 1 and 2 of the paper: the tunnel in 3D and from the top, and its cross-section with the detector layers.

    python -m curves.plot.tunnel_figures <output dir>

Drawn from grendel.geometry.grendel_geometry (the UP542 mirror), so they follow the geometry the scans use. The 3D
and top views are those of llpatcolliders higgs/tunnel_figure_panels.py (fiducial mesh, centreline, CMS IP); the
cross-section shows the tracker (two layers 24 cm apart on the roof and the left wall), the veto (floor and right
wall) and the fiducial inset, with an illustrative LLP decay.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from grendel.geometry import grendel_geometry as G

TRACKER, VETO, DAUGHTER = "#0072B2", "#D55E00", "#CC79A7"


def _equal_3d(ax, expand=(1.0, 1.0, 1.0)):
    ranges = []
    for (lo, hi), f, setter in zip((ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()), expand,
                                   (ax.set_xlim3d, ax.set_ylim3d, ax.set_zlim3d)):
        c, h = 0.5 * (lo + hi), 0.5 * (hi - lo) * f
        setter(c - h, c + h)
        ranges.append(2 * h)
    ax.set_box_aspect(ranges)


def tunnel_3d(out):
    m, p = G.mesh_fiducial, G.path_3d_fiducial
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")
    rng = np.random.default_rng(0)
    for fi in rng.choice(len(m.faces), min(1000, len(m.faces)), replace=False):
        t = m.vertices[m.faces[fi]]
        t = np.vstack([t, t[0]])
        ax.plot(t[:, 0], t[:, 2], t[:, 1], "b-", alpha=0.1, lw=0.5)
    ax.plot(p[:, 0], p[:, 2], p[:, 1], "r-", lw=3, label="Centreline")
    ax.scatter(0, 0, 0, color="green", s=200, label="CMS IP")
    a = 10
    for v, c, lab in (((a, 0, 0), "red", "X"), ((0, 0, a), "green", "Y (up)"), ((0, a, 0), "blue", "Z (beam)")):
        ax.quiver(0, 0, 0, *v, color=c, arrow_length_ratio=0.1)
        ax.text(*v, lab, fontsize=11)
    ax.set_xlabel("X (m)"); ax.set_ylabel("Z (beam, m)"); ax.set_zlabel("Y (up, m)")
    _equal_3d(ax, expand=(1.0, 1.2, 1.0))
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def tunnel_top(out):
    p = G.path_3d_fiducial
    half = G.tunnel_profile_points(inset=0)[:, 0].max()
    xz = p[:, [0, 2]]
    t = np.gradient(xz, axis=0)
    t /= np.linalg.norm(t, axis=1, keepdims=True)
    n = np.column_stack([-t[:, 1], t[:, 0]])
    band = np.vstack([xz + half * n, (xz - half * n)[::-1]])
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.fill(band[:, 0], band[:, 1], color="blue", alpha=0.3, label="Tunnel width")
    ax.plot(xz[:, 0], xz[:, 1], "r-", lw=2, label="Centreline")
    ax.scatter(0, 0, color="green", s=200, label="CMS IP")
    ax.set_xlabel("X (m)"); ax.set_ylabel("Z (beam, m)"); ax.set_title("Top View (X-Z Plane)")
    ax.grid(True, alpha=0.3); ax.axis("equal"); ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def cross_section(out):
    d = G.DETECTOR_THICKNESS
    nw, na = 16, 128
    outer, inner = (G.tunnel_profile_points(n_arch=na, n_wall=nw, inset=s) for s in (0.0, d))
    # vertex order: floor-left, floor-right, right wall up (nw), arch (na - 1), left wall down (nw)
    i_rtop, i_ltop = 1 + nw, 1 + nw + na          # right and left springline
    floor_y = outer[0, 1]
    fig, ax = plt.subplots(figsize=(8.5, 5.4))
    trk_o = np.vstack([outer[i_rtop:], outer[:1]])                    # arch + left wall, outer layer
    trk_i = np.vstack([inner[i_rtop:], inner[:1]])                    # inner layer, 24 cm in
    ax.plot(trk_o[:, 0], trk_o[:, 1], color=TRACKER, lw=4, label="Tracker surface")
    ax.plot(trk_i[:, 0], trk_i[:, 1], color=TRACKER, lw=4)
    veto = outer[:i_rtop + 1]                                         # floor + right wall
    ax.plot(veto[:, 0], veto[:, 1], color=VETO, lw=4, label="Scintillator veto")
    ax.plot([inner[1, 0], inner[i_rtop, 0]], [floor_y, inner[i_rtop, 1]], "--", color="0.6", lw=1.2,
            label=f"Fiducial inset ({d * 100:.0f} cm)")
    fid = np.vstack([inner[i_rtop:], inner[:1], [inner[1, 0], floor_y]])
    ax.fill(np.r_[fid[:, 0], inner[i_rtop, 0]], np.r_[fid[:, 1], inner[i_rtop, 1]], color="0.96", zorder=0)
    v = np.array([0.85, floor_y + 0.62])                               # illustrative decay
    ax.plot([v[0], v[0] + 1.6], [v[1], v[1] - 1.45], ":", color="0.2", lw=2.5, label="LLP (invisible)")
    for end in ((-0.33, outer[i_rtop + na // 2 - 22, 1]), (-1.45, 0.10)):
        ax.plot([v[0], end[0]], [v[1], end[1]], color=DAUGHTER, lw=2.2)
    ax.plot([], [], color=DAUGHTER, lw=2.2, label="SM daughters")
    ax.plot(*v, "o", color="0.15", ms=7)
    ax.text(0.30, 0.55, "SM", color=DAUGHTER, fontsize=12, weight="bold")
    ax.text(-0.60, -0.55, "SM", color=DAUGHTER, fontsize=12, weight="bold")
    ax.text(1.72, floor_y - 0.05, "LLP", rotation=-42, fontsize=12, weight="bold", color="0.2")
    ax.plot(0, 0, "+", color="k", ms=14, mew=1.5)
    ax.text(0, outer[:, 1].max() + 0.12, "Arch / Ceiling (tracker)", color=TRACKER, ha="center", fontsize=11, weight="bold")
    ax.text(-1.75, floor_y + 0.85, "Left wall (tracker)", color=TRACKER, rotation=90, ha="center", va="center",
            fontsize=11, weight="bold")
    ax.text(-0.05, floor_y + 0.36, "Floor (scintillator / veto)", color=VETO, ha="center", fontsize=10, weight="bold")
    ax.text(1.65, floor_y + 0.55, "Right wall\n(scintillator / veto)", color=VETO, fontsize=11, weight="bold")
    ax.text(-3.0, outer[:, 1].max() + 0.95, f"Two tracker layers\nseparated by {d * 100:.0f} cm", color=TRACKER,
            fontsize=13, weight="bold", va="top")
    gray = dict(arrowstyle="<->", color="0.5", lw=1)
    w = 2 * outer[i_rtop, 0]
    ax.annotate("", (-outer[i_rtop, 0], floor_y - 0.28), (outer[i_rtop, 0], floor_y - 0.28), arrowprops=gray)
    ax.text(0, floor_y - 0.50, f"{w:.2f} m (springline)", color="0.5", ha="center", fontsize=10)
    h = outer[:, 1].max() - floor_y
    ax.annotate("", (-2.05, floor_y), (-2.05, outer[:, 1].max()), arrowprops=gray)
    ax.text(-2.2, floor_y + h / 2, f"{h:.2f} m", color="0.5", rotation=90, ha="center", va="center", fontsize=10)
    ax.annotate("", (-2.7, floor_y), (-2.7, outer[i_ltop, 1]), arrowprops=gray)
    ax.text(-2.85, floor_y + (outer[i_ltop, 1] - floor_y) / 2, f"wall h = {outer[i_ltop, 1] - floor_y:.2f} m",
            color="0.5", rotation=90, ha="center", va="center", fontsize=10)
    ax.set_xlim(-3.15, 3.55); ax.set_ylim(floor_y - 0.65, outer[:, 1].max() + 1.0)
    ax.set_aspect("equal"); ax.grid(True, alpha=0.25)
    ax.set_xlabel("local x (m)", fontsize=13); ax.set_ylabel("local y (m)", fontsize=13)
    ax.legend(loc="upper right", fontsize=10.5, framealpha=0.95)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    o = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    o.mkdir(parents=True, exist_ok=True)
    tunnel_3d(o / "tunnel_3d.png")
    tunnel_top(o / "tunnel_topview.png")
    cross_section(o / "horseshoe_diagram.png")
    print(f"written to {o}: tunnel_3d.png tunnel_topview.png horseshoe_diagram.png")
