"""BC7 (HNL, muon dominance) from the PBC 2025 figure, raster route.

An independent digitisation of the same contours as ``vector_bc7``: the
PDF is rendered at 400 dpi, each contour is isolated by colour, thin
straight residuals (tick marks, leaders) are removed by a morphological
opening and a per-component PCA test, and the island edges are read off
column by column as the top and bottom of the coloured pixels.

    python -m curves.tools.raster_bc7
"""
from __future__ import annotations

import sys

import numpy as np

from ..plot.data import HNL_DATA
from .sources import BC7_PDF, require

OUT_DIR = HNL_DATA / "raster"
SOURCE = "Figures/5-PhysicsReach/BSM_benchmarks/pbc_bc7_gray_final.pdf"
FIG_LABEL = "fig:pbc_bc7_gray_final"
DATE = "2026-06-06"

# Pixel calibration of the plot frame at 400 dpi.
DPI = 400
LEFT, RIGHT, TOP, BOT = 442, 2922, 150, 1998
LOGO = (150, 430, 460, 830)   # inset logo, masked out

COLORS = {"ANUBIS": (75, 0, 130), "CODEX-b": (0, 139, 139), "FLArE": (231, 9, 0), "FASER2": (0, 100, 0),
          "SHiP": (218, 112, 214)}
# Tolerances chosen per curve to capture antialiasing.
TOL = {"ANUBIS": 70, "CODEX-b": 50, "FLArE": 80, "FASER2": 50, "SHiP": 60}
# Curve labels use the same hue; their pixel boxes are masked out.
LABEL_REGIONS = {"ANUBIS": (1130, 1330, 1750, 2700), "CODEX-b": (1280, 1450, 1700, 2330),
                 "FLArE": (1320, 1640, 1530, 2150), "FASER2": (1380, 1760, 1450, 2200),
                 "SHiP": (1600, 1990, 600, 1300)}
# Per-curve minimum mass to emit (what is visible in the figure; ANUBIS's
# dotted upper segment starts near 0.40 GeV).
M_START = {"ANUBIS": 0.40, "CODEX-b": 0.10, "FLArE": 0.10, "FASER2": 0.10, "SHiP": 0.10}


def render():
    import fitz
    from PIL import Image
    with fitz.open(str(require(BC7_PDF))) as doc:
        pix = doc[0].get_pixmap(dpi=DPI, alpha=False)
    return np.array(Image.frombytes("RGB", (pix.width, pix.height), pix.samples).convert("RGB"))


def px_to_data(px, py):
    lx = -1.0 + (np.asarray(px) - LEFT) / (RIGHT - LEFT) * 3.0
    ly = -2.0 + (np.asarray(py) - TOP) / (BOT - TOP) * (-10.0)
    return 10 ** lx, 10 ** ly


def color_mask(arr, name):
    target = np.array(COLORS[name], dtype=int)
    diff = np.abs(arr.astype(int) - target).sum(axis=2)
    m = diff <= TOL[name]
    y0, y1, x0, x1 = LOGO
    m[y0:y1, x0:x1] = False
    y0, y1, x0, x1 = LABEL_REGIONS[name]
    m[y0:y1, x0:x1] = False
    out = np.zeros_like(m)
    out[TOP:BOT + 1, LEFT:RIGHT + 1] = m[TOP:BOT + 1, LEFT:RIGHT + 1]
    return out


def is_thin_line(coords, min_n=20, min_len=25, max_ratio=0.012):
    if len(coords) < min_n:
        return False
    c = coords - coords.mean(axis=0)
    eig = np.sort(np.linalg.eigvalsh(np.cov(c.T)))[::-1]
    if eig[0] < 50:
        return False
    return eig[1] / eig[0] < max_ratio and np.sqrt(eig[0]) > min_len


def cleaned_mask(arr, name, open_size=3):
    from scipy.ndimage import binary_opening, label
    m1 = binary_opening(color_mask(arr, name), structure=np.ones((open_size, open_size)))
    lab, n = label(m1, structure=np.ones((3, 3)))
    out = m1.copy()
    for i in range(1, n + 1):
        coords = np.argwhere(lab == i)
        if is_thin_line(coords):
            out[coords[:, 0], coords[:, 1]] = False
    return out


def envelope(mask, x_window=15, gap_thresh=30):
    """Per column: (top, bottom, n_clusters) of the coloured pixels pooled
    over a small x window, split into clusters at vertical gaps."""
    ys, xs = np.where(mask)
    by_x: dict[int, list[int]] = {}
    for x, y in zip(xs, ys):
        by_x.setdefault(int(x), []).append(int(y))
    env = {}
    for x in range(LEFT, RIGHT + 1):
        pooled = []
        for dx in range(-x_window, x_window + 1):
            if x + dx in by_x:
                pooled.extend(by_x[x + dx])
        if not pooled:
            continue
        pooled.sort()
        clusters = [[pooled[0]]]
        for y in pooled[1:]:
            if y - clusters[-1][-1] > gap_thresh:
                clusters.append([y])
            else:
                clusters[-1].append(y)
        centres = sorted(float(np.mean(c)) for c in clusters)
        env[x] = (centres[0], centres[-1], len(centres))
    return env


def build_mass_grid():
    # ~830 px/decade, so a sub-percent mass step is ~2 px: the pooling
    # window smooths column noise without losing real structure.
    base = np.logspace(np.log10(0.1), np.log10(30), 300)
    dense = np.concatenate([np.arange(0.25, 0.45 + 1e-9, 0.002), np.arange(1.60, 2.00 + 1e-9, 0.002),
                            np.arange(4.80, 5.40 + 1e-9, 0.005)])
    return np.sort(np.unique(np.round(np.concatenate([base, dense]), 5)))


MASS_GRID = build_mass_grid()


def query(env, x_target, x_window=20, mode="median"):
    yts, ybs, ncs = [], [], []
    for dx in range(-x_window, x_window + 1):
        if x_target + dx in env:
            yt, yb, n = env[x_target + dx]
            yts.append(yt)
            ybs.append(yb)
            ncs.append(n)
    if not yts:
        return None
    if mode == "extrema":
        return float(min(yts)), float(max(ybs)), int(max(ncs))
    return float(np.median(yts)), float(np.median(ybs)), int(np.median(ncs))


def extract_rows(mask, env, name, query_mode="median"):
    ys, xs = np.where(mask)
    if not len(xs):
        return []
    xmin, xmax = int(xs.min()), int(xs.max())
    rows = []
    for m in MASS_GRID:
        if m < M_START[name] - 1e-6:
            continue
        x_target = int(round(LEFT + (np.log10(m) + 1.0) / 3.0 * (RIGHT - LEFT)))
        if x_target < xmin - 30 or x_target > xmax + 30:
            continue
        q = query(env, x_target, x_window=18, mode=query_mode)
        if q is None:
            continue
        y_top, y_bot, _ = q
        _, u_upper = px_to_data(0, y_top)
        _, u_lower = px_to_data(0, y_bot)
        if u_lower > u_upper:
            u_lower, u_upper = u_upper, u_lower
        rows.append((float(m), float(u_lower), float(u_upper)))
    return rows


def repair_ship_upper(rows):
    """Fill sparse-raster gaps in the SHiP high-U branch by log
    interpolation from neighbouring visible upper-edge pixels, instead of
    emitting artificial sentinel jumps."""
    arr = np.array(rows, dtype=float)
    if len(arr) < 2:
        return rows
    good = (arr[:, 2] > 1e-5) & (arr[:, 2] < 1e-2) & (arr[:, 0] < 3.7)
    if good.sum() < 2:
        return rows
    lx = np.log10(arr[good, 0])
    ly = np.log10(arr[good, 2])
    order = np.argsort(lx)
    lx, ly = lx[order], ly[order]
    missing = (arr[:, 2] < 1e-5) & (arr[:, 0] < 3.7)
    before = missing & (arr[:, 0] < 10 ** lx.min())
    inside = missing & (arr[:, 0] >= 10 ** lx.min()) & (arr[:, 0] <= 10 ** lx.max())
    arr[before, 2] = 10 ** ly[0]
    arr[inside, 2] = 10 ** np.interp(np.log10(arr[inside, 0]), lx, ly)
    return [tuple(row) for row in arr]


def write_dat(name, fname, rows, extra=None):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / fname
    with path.open("w") as fh:
        fh.write(f"# {name} HNL |Umu|^2 sensitivity (BC7 benchmark)\n")
        fh.write(f"# Digitized from arXiv:2505.00947v2, {FIG_LABEL}, {DATE}\n")
        fh.write(f"# Source file: {SOURCE}\n")
        for ln in extra or ():
            fh.write(f"# {ln}\n")
        fh.write("# mass_GeV  u2_min  u2_max\n")
        for m, lo, hi in rows:
            fh.write(f"{m:.4f}  {lo:.3e}  {hi:.3e}\n")
    print(f"wrote {path} ({len(rows)} rows)")


def main() -> int:
    arr = render()
    masks = {name: cleaned_mask(arr, name) for name in COLORS}
    envs = {name: envelope(mask) for name, mask in masks.items()}
    for name, mask in masks.items():
        ys, xs = np.where(mask)
        if len(xs):
            print(f"  {name}: {mask.sum()} px, m={px_to_data(xs.min(), 0)[0]:.3g}..{px_to_data(xs.max(), 0)[0]:.3g}")
    out = {name: extract_rows(masks[name], envs[name], name) for name in ("ANUBIS", "CODEX-b", "FLArE", "FASER2")}
    out["SHiP"] = repair_ship_upper(extract_rows(masks["SHiP"], envs["SHiP"], "SHiP", query_mode="extrema"))
    write_dat("ANUBIS", "ANUBIS_Umu.dat", out["ANUBIS"])
    write_dat("CODEX-b", "CODEX-b_Umu.dat", out["CODEX-b"])
    write_dat("FLArE", "FLArE_Umu.dat", out["FLArE"],
              extra=["FLArE and FASER2 contours overlap to within line width across",
                     "most of the mass range; values are co-located with FASER2_Umu.dat."])
    write_dat("FASER2", "FASER2_Umu.dat", out["FASER2"],
              extra=["FASER2 and FLArE contours overlap to within line width across",
                     "most of the mass range; values are co-located with FLArE_Umu.dat."])
    write_dat("SHiP", "SHiP_Umu.dat", out["SHiP"],
              extra=["The SHiP upper branch is sampled with an extrema window and",
                     "short sparse-raster gaps are log-interpolated from neighboring",
                     "visible upper-edge pixels to avoid artificial branch toggles."])
    return 0


if __name__ == "__main__":
    sys.exit(main())
