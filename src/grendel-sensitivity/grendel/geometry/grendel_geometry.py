"""
GRENDEL detector geometry module.

Provides the fiducial volume mesh and ray-casting utilities for the
GRENDEL tunnel detector above CMS IP5.

Coordinate convention (CMS physicist frame):
    X = horizontal transverse (+ towards the LHC centre)
    Y = up (perpendicular to the beam)
    Z = beam axis

The tunnel is the UP542 gallery of the GRENDEL geometry study (grendel-geometry, since October
2026), taken unchanged from its single source, https://gitlab.cern.ch/grendel/grendel-geometry
(results/geometry/), mirrored in grendel/geometry/data/up542/ and checked against its
manifest on import (update with scripts/sync_up542_geometry.py; see GEOMETRY_UP542.md):
the centreline with the finished floor level per chainage from the CERN civil drawings
(re-issued "as built" by the contractor in 2004), the cross-section from the 2015 laser scan
by the CMS survey team (2.98 m wide, 2.95 m high, circular roof R 1.52 m), the
correct handedness and the true length (109.76 m), in the official CMS physicist frame
(the level survey frame pitched by 12.34 mrad = 0.707 deg; CMS-I-UR-0281).
Cross-section faces are perpendicular to the centreline, with "up" along the true vertical.

The tunnel of releases v1-v2.2 (the curves of arXiv:2609.00152v1/v2: 47-point centreline
at Y = 22 m, 2.90 m x 3.15 m half-ellipse) is in the git tags v1 ... v2.2.
"""

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import trimesh
from tqdm import tqdm

# ============================================================
# Physical constants
# ============================================================
SPEED_OF_LIGHT = 299792458.0  # m/s

# ============================================================
# The UP542 geometry: mirror of the single source (never edit it here)
# ============================================================
UP542_DIR = Path(__file__).resolve().parent / "data" / "up542"


def _check_up542_mirror():
    man = json.loads((UP542_DIR / "manifest.json").read_text())
    bad = [f for f, h in man["files"].items()
           if hashlib.sha256((UP542_DIR / f).read_bytes()).hexdigest() != h]
    if bad:
        raise RuntimeError(f"grendel/geometry/data/up542 differs from its manifest ({', '.join(bad)}): "
                           "re-sync it with scripts/sync_up542_geometry.py instead of editing it")
    return man["version"]


UP542_VERSION = _check_up542_mirror()
_FRAME = json.loads((UP542_DIR / "UP542_frame.json").read_text())
_SECTION = json.loads((UP542_DIR / "UP542_section.json").read_text())["measured"]

# Official CMS frame: the level survey frame pitched about the transverse axis (CMS-I-UR-0281,
# EDMS 2028417); the beam (Z_CMS) rises towards L54.
CMS_PITCH_RAD = _FRAME["cms_physicist_frame"]["pitch_mrad"] * 1e-3
# True vertical in CMS coordinates (the level frame's Z axis): (0, cos a, sin a)
WORLD_UP = np.array([0.0, np.cos(CMS_PITCH_RAD), np.sin(CMS_PITCH_RAD)])

# ============================================================
# Tunnel cross-section: from the 2015 laser scan of UP542 by the CMS survey team (CMS-I-UR-0261;
# cloud re-fitted in 2017, tilt removed), medians over the scanned sections: straight walls from the
# floor to the springline, circular roof (2.882 / 1.654 / 2.947 / 1.516 m)
# ============================================================
TUNNEL_FLOOR_WIDTH = _SECTION["floor_width_m"]          # at the floor, above the kerbs (m)
TUNNEL_SPRING_H = _SECTION["springline_above_floor_m"]  # springline above the floor (m)
TUNNEL_CROWN_H = _SECTION["crown_above_floor_m"]        # floor to crown (m)
TUNNEL_ROOF_R = _SECTION["roof_radius_m"]               # roof circle radius (m)
TUNNEL_CENTROID_H = _SECTION["centroid_above_floor_m"]  # section centroid above the floor (m)

# Single source of truth for the detector shell radial extent. This sets the
# fiducial-volume inset AND the tracker-layer spacing (inner layer at the
# fiducial face, outer layer at the tunnel wall). Change it in this ONE place
# to study a different detector thickness; the mesh and all reconstruction
# (signal and background) pick it up consistently.
DETECTOR_THICKNESS = 0.24  # 24 cm


# ============================================================
# Tunnel profile (2D cross-section)
# ============================================================

def tunnel_profile_points(n_arch=32, n_wall=4, inset=0.0, inset_floor=False):
    """
    Generate 2D cross-section points for the tunnel profile.

    Parameters
    ----------
    n_arch : int
        Number of points along the arch.
    n_wall : int
        Number of points along each wall segment.
    inset : float
        Inward offset (m) from the tunnel walls (e.g. detector thickness).
    inset_floor : bool
        Whether to also inset the floor.

    Returns
    -------
    points : ndarray, shape (N, 2)
        2D profile coordinates centred on the cross-section centroid.
        Vertex order: floor-left, floor-right, right wall up (n_wall),
        arch (n_arch - 1), left wall down (n_wall).
    """
    R, yc, ys = TUNNEL_ROOF_R, TUNNEL_CROWN_H - TUNNEL_ROOF_R, TUNNEL_SPRING_H
    half_spring = np.sqrt(R ** 2 - (ys - yc) ** 2) - inset
    half_floor = TUNNEL_FLOOR_WIDTH / 2 - inset
    floor_y = inset if inset_floor else 0.0
    wall = [[half_floor + (half_spring - half_floor) * i / n_wall, floor_y + (ys - floor_y) * i / n_wall]
            for i in range(1, n_wall + 1)]
    phi0 = np.arctan2(ys - yc, half_spring)                 # top of the (inset) wall on the roof circle
    arch = [[(R - inset) * np.cos(phi), yc + (R - inset) * np.sin(phi)]
            for phi in phi0 + (np.pi - 2 * phi0) * np.arange(1, n_arch) / n_arch]
    points = np.array([[-half_floor, floor_y], [half_floor, floor_y]] + wall + arch
                      + [[-x, y] for x, y in reversed(wall)], dtype=float)
    # centred on the centroid of the full (not inset) section: the source's value, so that the
    # profile sits exactly on the source's centreline
    points[:, 1] -= TUNNEL_CENTROID_H
    return points


# ============================================================
# Mesh extrusion along a 3D path
# ============================================================

def create_profile_mesh(path_points, profile_2d):
    """
    Extrude a 2D profile along a 3D path to create a triangulated mesh.

    The world-up reference vector is WORLD_UP: the true vertical in the CMS physicist frame (pitched 12.34 mrad).
    """
    n_profile = len(profile_2d)
    vertices = []
    faces = []

    for i in range(len(path_points)):
        if i == 0:
            tangent = path_points[1] - path_points[0]
        elif i == len(path_points) - 1:
            tangent = path_points[i] - path_points[i - 1]
        else:
            tangent = path_points[i + 1] - path_points[i - 1]
        tangent = tangent / np.linalg.norm(tangent)

        # up = the true vertical (CMS Y, pitched for the grendel-geometry tunnel)
        if abs(tangent[1]) < 0.9:
            world_up = WORLD_UP
        else:
            world_up = np.array([0, 0, 1])

        right = np.cross(tangent, world_up)
        right = right / np.linalg.norm(right)
        up = np.cross(right, tangent)
        up = up / np.linalg.norm(up)

        for j in range(n_profile):
            offset = profile_2d[j, 0] * right + profile_2d[j, 1] * up
            vertices.append(path_points[i] + offset)

        if i > 0:
            for j in range(n_profile):
                v1 = (i - 1) * n_profile + j
                v2 = (i - 1) * n_profile + (j + 1) % n_profile
                v3 = i * n_profile + (j + 1) % n_profile
                v4 = i * n_profile + j
                faces.append([v1, v4, v3])
                faces.append([v1, v3, v2])

    # Cap the start
    center_start = len(vertices)
    vertices.append(path_points[0].copy())
    for j in range(n_profile):
        faces.append([center_start, (j + 1) % n_profile, j])

    # Cap the end
    center_end = len(vertices)
    vertices.append(path_points[-1].copy())
    last = (len(path_points) - 1) * n_profile
    for j in range(n_profile):
        faces.append([center_end, last + j, last + (j + 1) % n_profile])

    return np.array(vertices), np.array(faces)


# ============================================================
# Direction helpers
# ============================================================

def eta_phi_to_direction(eta, phi):
    """Convert (η, φ) to a unit direction vector in CMS coordinates."""
    theta = 2 * np.arctan(np.exp(-eta))
    dx = np.sin(theta) * np.cos(phi)
    dy = np.sin(theta) * np.sin(phi)
    dz = np.cos(theta)
    return np.array([dx, dy, dz])


def calculate_decay_length(momentum, mass, lifetime):
    """Compute lab-frame decay length βγcτ (metres)."""
    energy = np.sqrt(momentum**2 + mass**2)
    beta = momentum / energy
    gamma = energy / mass
    return gamma * beta * SPEED_OF_LIGHT * lifetime


# ============================================================
# Geometry caching (ray-cast each particle once)
# ============================================================

def cache_geometry(csv_file, mesh, origin):
    """
    Ray-cast every particle in *csv_file* against *mesh* and return
    a dict of cached entry/exit distances and kinematic quantities.
    """
    df = pd.read_csv(csv_file)
    n = len(df)
    origin = np.array(origin)

    hits = np.zeros(n, dtype=bool)
    entry_d = np.full(n, np.nan)
    exit_d = np.full(n, np.nan)
    direction_arr = np.zeros((n, 3))
    momentum = df['momentum'].values
    mass = df['mass'].values
    energy = np.sqrt(momentum**2 + mass**2)
    gamma = energy / mass
    beta = momentum / energy

    print(f"Caching fiducial volume geometry for {n} particles...")
    for idx, row in tqdm(df.iterrows(), total=n, desc="Ray-casting"):
        direction = eta_phi_to_direction(row['eta'], row['phi'])
        direction_arr[idx] = direction
        locations, _, _ = mesh.ray.intersects_location(
            ray_origins=[origin], ray_directions=[direction])
        if len(locations) >= 2:
            hits[idx] = True
            distances = sorted(
                [np.linalg.norm(loc - origin) for loc in locations])
            entry_d[idx] = distances[0]
            exit_d[idx] = distances[1]

    n_hits = hits.sum()
    print(f"  {n_hits} / {n} particles hit fiducial volume "
          f"({n_hits / n * 100:.1f}%)")
    if n_hits > 0:
        print(f"  Mean path length: "
              f"{(exit_d[hits] - entry_d[hits]).mean():.2f} m")

    return {
        'hits': hits, 'entry_d': entry_d, 'exit_d': exit_d,
        'gamma': gamma, 'beta': beta, 'momentum': momentum, 'mass': mass,
        'direction': direction_arr, 'event': df['event'].values,
    }


# ============================================================
# Tunnel centreline (CERN civil drawings, 2004; from data/up542/UP542_centreline.csv)
# ============================================================
TUNNEL_CENTRELINE_CSV = UP542_DIR / "UP542_centreline.csv"
CENTRELINE_STEP = 4       # use every 4th row of the CSV (~1 m; chord error < 2 cm on the R = 8 m bends)


def tunnel_centreline(y_position=None, step=CENTRELINE_STEP):
    """
    Section centroids of the gallery in the CMS physicist frame (m), and their chainage.
    y_position=None: the source's centroid (finished floor level per chainage + the scanned section's
    centroid height, pitched into the CMS frame); a number: a flat centroid at that height above IP5
    (level frame, then pitched), for comparisons.
    """
    with open(TUNNEL_CENTRELINE_CSV, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(line for line in fh if not line.startswith("#")))
    idx = np.unique(np.r_[np.arange(0, len(rows), step), len(rows) - 1])
    col = lambda k: np.array([float(rows[i][k]) for i in idx])
    if y_position is None:
        return np.c_[col("x_CMS_m"), col("y_CMS_m"), col("z_CMS_m")], col("chainage_m")
    xb, yt = col("Xb_m"), col("Yt_m")
    zc = np.full_like(xb, float(y_position))
    a = CMS_PITCH_RAD
    return np.c_[-yt, zc * np.cos(a) + xb * np.sin(a), -xb * np.cos(a) + zc * np.sin(a)], col("chainage_m")


# Mean height of the tunnel centroid above IP5 (informational; the floor varies along the gallery)
Y_POSITION = float(np.mean(tunnel_centreline()[0][:, 1]))


# ============================================================
# Build the default fiducial-volume mesh
# ============================================================

def build_fiducial_mesh(y_position=None,
                        detector_thickness=DETECTOR_THICKNESS):
    """
    Construct the fiducial-volume trimesh for the GRENDEL tunnel.

    Parameters
    ----------
    y_position : float or None
        None = floor level per chainage (default); a number = flat centroid height above IP5 (m).
    detector_thickness : float
        Inset from the tunnel wall defining the fiducial volume (m).

    Returns
    -------
    mesh : trimesh.Trimesh
        Closed triangulated mesh of the fiducial volume.
    path_3d : ndarray, shape (N, 3)
        3D tunnel centreline coordinates.
    """
    path_3d, _ = tunnel_centreline(y_position)
    # Walked from the south end towards PX56. The wall roles (TRACKER_SURFACES,
    # SCINTILLATOR_SURFACES) are left/right of the walking direction; the old tunnel
    # was the mirror image walked from PX56, with the veto ("Right Wall") facing IP5
    # on 80 % of the gallery. Walking the true gallery the other way keeps that.
    path_3d = path_3d[::-1]
    profile = tunnel_profile_points(inset=detector_thickness,
                                    inset_floor=False)
    verts, faces = create_profile_mesh(path_3d, profile)
    mesh = trimesh.Trimesh(vertices=verts, faces=faces)
    # The Frenet-frame extrusion in create_profile_mesh produces inconsistent
    # face winding along curved sections; without fix_normals, mesh.volume
    # collapses to roughly 1/3 of the true volume (signed integral cancels).
    mesh.fix_normals()

    print(f"Fiducial volume built (UP542, grendel-geometry version {UP542_VERSION}):")
    print(f"  Volume:    {mesh.volume:.1f} m³")
    print(f"  Y offset:  {'floor per chainage, mean %.2f' % path_3d[:, 1].mean() if y_position is None else y_position} m")
    print(f"  Thickness: {detector_thickness * 100:.0f} cm")

    return mesh, path_3d


# Convenience: build the default mesh on import
mesh_fiducial, path_3d_fiducial = build_fiducial_mesh()


# ============================================================
# Surface classification (shared with signal_surface_hitmap.py)
# ============================================================
# The tunnel wall is divided into four named surfaces. In the GRENDEL
# design the Floor and Right Wall are instrumented with scintillator
# (hit/timing only, no tracks) while the Arch/Ceiling and Left Wall carry
# the tracking layers. A two-body decay is only reconstructable if BOTH
# daughters land on a tracker surface.

# Centreline arc-length parametrisation
seg_lengths = np.array(
    [np.linalg.norm(path_3d_fiducial[i + 1] - path_3d_fiducial[i])
     for i in range(len(path_3d_fiducial) - 1)])
cumulative_length = np.concatenate([[0], np.cumsum(seg_lengths)])
total_length = cumulative_length[-1]

# Profile-angle boundaries between the four surfaces, derived from the
# actual profile vertices (default n_arch=32, n_wall=4).
_profile_pts = tunnel_profile_points(inset=0.0)
_profile_angles = np.arctan2(_profile_pts[:, 1], _profile_pts[:, 0])
_profile_angles = np.where(_profile_angles < 0,
                           _profile_angles + 2 * np.pi, _profile_angles)

# Vertex order from tunnel_profile_points:
#   [floor-left, floor-right, right wall up..., arch..., left wall down...]
_N_WALL = 4
_N_ARCH_PTS = 31  # n_arch - 1 interior points for n_arch=32
theta_floor_left  = _profile_angles[0]
theta_floor_right = _profile_angles[1]
theta_rwall_top   = _profile_angles[1 + _N_WALL]
theta_arch_top    = _profile_angles[1 + _N_WALL + (_N_ARCH_PTS // 2)]
theta_lwall_top   = _profile_angles[1 + _N_WALL + _N_ARCH_PTS]

TRACKER_SURFACES = ('Arch/Ceiling', 'Left Wall')
SCINTILLATOR_SURFACES = ('Floor', 'Right Wall')


# ============================================================
# Profile arc-length parametrization (around the cross-section)
# ============================================================
# Cumulative arc-length along the closed polygon walk of the profile
# vertices. The walk is CCW so the unwrapped angle increases by 2π
# over one full traversal.
_seg_lens = np.linalg.norm(
    np.diff(_profile_pts, axis=0, append=_profile_pts[:1]), axis=1)
_profile_cum_s = np.concatenate([[0], np.cumsum(_seg_lens)])  # length N+1
profile_perimeter = float(_profile_cum_s[-1])

# Build a monotonically-increasing (theta, s) lookup by unwrapping the
# polygon angles. Close the loop by repeating the first vertex at
# theta + 2π.
_theta_walk = _profile_angles.copy()
for i in range(1, len(_theta_walk)):
    while _theta_walk[i] < _theta_walk[i - 1]:
        _theta_walk[i] += 2 * np.pi
_theta_walk = np.concatenate([_theta_walk, [_theta_walk[0] + 2 * np.pi]])
_s_walk = np.concatenate([_profile_cum_s[:-1], [profile_perimeter]])


def arc_length_at_theta(theta):
    """
    Cumulative arc-length s(theta) along the cross-section profile (m),
    measured from the start of the polygon walk (floor-left corner).
    Accepts scalar or array inputs in [0, 2π).
    """
    theta = np.atleast_1d(np.asarray(theta, dtype=float))
    theta_mod = theta % (2 * np.pi)
    base = _theta_walk[0]
    theta_shifted = np.where(theta_mod < base, theta_mod + 2 * np.pi, theta_mod)
    return np.interp(theta_shifted, _theta_walk, _s_walk)


def arc_distance_between_thetas(theta_a, theta_b):
    """
    Shortest arc-length distance (m) along the profile contour between
    two points at profile angles theta_a and theta_b. Naturally bounded
    by profile_perimeter / 2.
    """
    s_a = arc_length_at_theta(theta_a)
    s_b = arc_length_at_theta(theta_b)
    diff = np.abs(s_a - s_b)
    return np.minimum(diff, profile_perimeter - diff)


def _in_arc(theta, a, b):
    """True where theta is in the half-open CCW arc [a, b) (handles wrap)."""
    theta = np.asarray(theta, dtype=float) % (2 * np.pi)
    if a <= b:
        return (theta >= a) & (theta < b)
    return (theta >= a) | (theta < b)


def classify_by_theta(theta):
    """Map a scalar profile angle to a surface name."""
    if bool(_in_arc(theta, theta_floor_right, theta_rwall_top)):
        return 'Right Wall'
    if bool(_in_arc(theta, theta_rwall_top, theta_lwall_top)):
        return 'Arch/Ceiling'
    if bool(_in_arc(theta, theta_lwall_top, theta_floor_left)):
        return 'Left Wall'
    return 'Floor'


def classify_exit_point(point, path_3d, cumulative_length):
    """
    For a 3D point on the tunnel wall, return (s, theta, x_local, y_local):
    arc-length s along the centreline, profile angle theta, and the local
    transverse coordinates at the nearest centreline segment.
    """
    best_s = 0.0
    best_dist_sq = np.inf
    best_x = 0.0
    best_y = 0.0

    for i in range(len(path_3d) - 1):
        seg = path_3d[i + 1] - path_3d[i]
        seg_len = np.linalg.norm(seg)
        if seg_len == 0:
            continue
        seg_hat = seg / seg_len
        t = np.clip(np.dot(point - path_3d[i], seg_hat), 0, seg_len)
        closest = path_3d[i] + t * seg_hat
        diff = point - closest
        dist_sq = np.dot(diff, diff)

        if dist_sq < best_dist_sq:
            best_dist_sq = dist_sq
            best_s = cumulative_length[i] + t
            tangent = seg_hat
            if abs(tangent[1]) < 0.9:
                world_up = WORLD_UP
            else:
                world_up = np.array([0., 0., 1.])
            right = np.cross(tangent, world_up)
            right /= np.linalg.norm(right)
            up = np.cross(right, tangent)
            up /= np.linalg.norm(up)
            best_x = np.dot(diff, right)
            best_y = np.dot(diff, up)

    theta = np.arctan2(best_y, best_x)
    if theta < 0:
        theta += 2 * np.pi

    return best_s, theta, best_x, best_y


def classify_points(points):
    """
    Vectorised surface mapping for many wall points at once.

    Parameters
    ----------
    points : ndarray, shape (M, 3)

    Returns
    -------
    theta : ndarray, shape (M,)  — profile angle of each point
    s     : ndarray, shape (M,)  — arc-length along the centreline
    """
    points = np.asarray(points, dtype=float)
    m = len(points)
    best_d2 = np.full(m, np.inf)
    best_x = np.zeros(m)
    best_y = np.zeros(m)
    best_s = np.zeros(m)

    for i in range(len(path_3d_fiducial) - 1):
        seg = path_3d_fiducial[i + 1] - path_3d_fiducial[i]
        seg_len = np.linalg.norm(seg)
        if seg_len == 0:
            continue
        seg_hat = seg / seg_len
        if abs(seg_hat[1]) < 0.9:
            world_up = WORLD_UP
        else:
            world_up = np.array([0., 0., 1.])
        right = np.cross(seg_hat, world_up)
        right /= np.linalg.norm(right)
        up = np.cross(right, seg_hat)
        up /= np.linalg.norm(up)

        rel = points - path_3d_fiducial[i]
        t = np.clip(rel @ seg_hat, 0, seg_len)
        closest = path_3d_fiducial[i] + np.outer(t, seg_hat)
        diff = points - closest
        d2 = np.einsum('ij,ij->i', diff, diff)
        upd = d2 < best_d2
        best_d2[upd] = d2[upd]
        best_x[upd] = diff[upd] @ right
        best_y[upd] = diff[upd] @ up
        best_s[upd] = cumulative_length[i] + t[upd]

    theta = np.arctan2(best_y, best_x)
    theta = np.where(theta < 0, theta + 2 * np.pi, theta)
    return theta, best_s


def classify_points_with_basis(points):
    """
    Same surface mapping as classify_points, but also returns the local
    cavern basis at each point: (tangent, right, up) where ``tangent``
    is the centreline tangent at the nearest segment ("down the tunnel"),
    and (right, up) span the profile plane. The radial outward normal at
    profile angle theta is cos(theta)*right + sin(theta)*up; the arc
    tangent (around the cross-section) is -sin(theta)*right + cos(theta)*up.

    Returns
    -------
    theta   : ndarray, shape (M,)
    s       : ndarray, shape (M,)
    tangent : ndarray, shape (M, 3)
    right   : ndarray, shape (M, 3)
    up      : ndarray, shape (M, 3)
    """
    points = np.asarray(points, dtype=float)
    m = len(points)
    best_d2 = np.full(m, np.inf)
    best_x = np.zeros(m)
    best_y = np.zeros(m)
    best_s = np.zeros(m)
    best_tangent = np.zeros((m, 3))
    best_right = np.zeros((m, 3))
    best_up = np.zeros((m, 3))

    for i in range(len(path_3d_fiducial) - 1):
        seg = path_3d_fiducial[i + 1] - path_3d_fiducial[i]
        seg_len = np.linalg.norm(seg)
        if seg_len == 0:
            continue
        seg_hat = seg / seg_len
        if abs(seg_hat[1]) < 0.9:
            world_up = WORLD_UP
        else:
            world_up = np.array([0., 0., 1.])
        right = np.cross(seg_hat, world_up)
        right /= np.linalg.norm(right)
        up = np.cross(right, seg_hat)
        up /= np.linalg.norm(up)

        rel = points - path_3d_fiducial[i]
        t = np.clip(rel @ seg_hat, 0, seg_len)
        closest = path_3d_fiducial[i] + np.outer(t, seg_hat)
        diff = points - closest
        d2 = np.einsum('ij,ij->i', diff, diff)
        upd = d2 < best_d2
        best_d2[upd] = d2[upd]
        best_x[upd] = diff[upd] @ right
        best_y[upd] = diff[upd] @ up
        best_s[upd] = cumulative_length[i] + t[upd]
        best_tangent[upd] = seg_hat
        best_right[upd] = right
        best_up[upd] = up

    theta = np.arctan2(best_y, best_x)
    theta = np.where(theta < 0, theta + 2 * np.pi, theta)
    return theta, best_s, best_tangent, best_right, best_up


def local_transverse_xy(points):
    """
    Local cross-section coordinates (x_local, y_local) of each (M,3) point:
    the offset from the nearest centreline point projected onto (right, up).
    Same projection as classify_points, but returns the 2-D profile coords.
    """
    pts = np.asarray(points, dtype=float)
    m = len(pts)
    best_d2 = np.full(m, np.inf)
    bx = np.zeros(m)
    by = np.zeros(m)
    for i in range(len(path_3d_fiducial) - 1):
        seg = path_3d_fiducial[i + 1] - path_3d_fiducial[i]
        seg_len = np.linalg.norm(seg)
        if seg_len == 0:
            continue
        seg_hat = seg / seg_len
        wu = WORLD_UP if abs(seg_hat[1]) < 0.9 else np.array([0., 0., 1.])
        right = np.cross(seg_hat, wu); right /= np.linalg.norm(right)
        up = np.cross(right, seg_hat); up /= np.linalg.norm(up)
        rel = pts - path_3d_fiducial[i]
        t = np.clip(rel @ seg_hat, 0, seg_len)
        diff = pts - (path_3d_fiducial[i] + np.outer(t, seg_hat))
        d2 = np.einsum('ij,ij->i', diff, diff)
        upd = d2 < best_d2
        best_d2[upd] = d2[upd]
        bx[upd] = diff[upd] @ right
        by[upd] = diff[upd] @ up
    return bx, by


# Fast point-in-fiducial test: project to the nearest centreline and check the
# local (x,y) against the 2-D fiducial profile polygon. ~180x faster than
# trimesh mesh.contains (which ray-casts per point) and consistent with the
# rest of the centreline-based geometry.
from matplotlib.path import Path as _MplPath
_fid_profile_2d = tunnel_profile_points(inset=DETECTOR_THICKNESS)
_fid_path_2d = _MplPath(np.vstack([_fid_profile_2d, _fid_profile_2d[:1]]))


def points_in_fiducial(points):
    """Bool array: True where each (M,3) point lies inside the fiducial volume.
    Projects to the nearest centreline segment, tests the local cross-section
    against the inset profile, and rejects points past the end-caps (nearest
    projection clamped at the first/last node). Fast, exact proxy for
    mesh.contains away from machine-precision boundaries."""
    pts = np.asarray(points, dtype=float)
    m = len(pts)
    best_d2 = np.full(m, np.inf)
    bx = np.zeros(m)
    by = np.zeros(m)
    beyond = np.zeros(m, dtype=bool)
    nseg = len(path_3d_fiducial) - 1
    for i in range(nseg):
        seg = path_3d_fiducial[i + 1] - path_3d_fiducial[i]
        seg_len = np.linalg.norm(seg)
        if seg_len == 0:
            continue
        seg_hat = seg / seg_len
        wu = WORLD_UP if abs(seg_hat[1]) < 0.9 else np.array([0., 0., 1.])
        right = np.cross(seg_hat, wu); right /= np.linalg.norm(right)
        up = np.cross(right, seg_hat); up /= np.linalg.norm(up)
        rel = pts - path_3d_fiducial[i]
        t_raw = rel @ seg_hat
        t = np.clip(t_raw, 0, seg_len)
        diff = pts - (path_3d_fiducial[i] + np.outer(t, seg_hat))
        d2 = np.einsum('ij,ij->i', diff, diff)
        upd = d2 < best_d2
        best_d2[upd] = d2[upd]
        bx[upd] = diff[upd] @ right
        by[upd] = diff[upd] @ up
        # past the end-caps: nearest projection clamped at the global ends
        if i == 0:
            beyond[upd] = t_raw[upd] < 0
        elif i == nseg - 1:
            beyond[upd] = t_raw[upd] > seg_len
        else:
            beyond[upd] = False
    in_profile = _fid_path_2d.contains_points(np.column_stack([bx, by]))
    return in_profile & ~beyond


def points_on_tracker(points):
    """
    Bool array: True where each (M, 3) wall point lands on a tracker
    surface (Arch/Ceiling or Left Wall). The tracker spans the single
    CCW arc [theta_rwall_top, theta_floor_left); the complement
    (Floor + Right Wall) is scintillator.
    """
    theta, _ = classify_points(points)
    return _in_arc(theta, theta_rwall_top, theta_floor_left)
