"""
GRENDEL detector geometry module.

Provides the fiducial volume mesh and ray-casting utilities for the
GRENDEL tunnel detector above CMS IP5.

Coordinate convention (CMS standard):
    X = horizontal transverse
    Y = vertical (up)
    Z = beam axis

The tunnel centerline lies in the X-Z plane at fixed Y = Y_POSITION.
Cross-section faces are oriented in the transverse X-Y plane.
"""

import numpy as np
import pandas as pd
import trimesh
from tqdm import tqdm

# ============================================================
# Physical constants
# ============================================================
SPEED_OF_LIGHT = 299792458.0  # m/s

# ============================================================
# Tunnel cross-section definition (Position C measurements)
# ============================================================
TUNNEL_ALPHA = 2.90   # floor width (m)
TUNNEL_BETA  = 3.15   # total height (m)
TUNNEL_GAMMA = 2.90   # arch width at springline (m)
TUNNEL_DELTA = 1.90   # arch height (m)
TUNNEL_WALL_HEIGHT = TUNNEL_BETA - TUNNEL_DELTA  # 1.25 m

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
    """
    half_w = TUNNEL_GAMMA / 2 - inset
    half_floor = TUNNEL_ALPHA / 2 - inset
    wall_h = TUNNEL_WALL_HEIGHT
    a = TUNNEL_GAMMA / 2 - inset
    b = TUNNEL_DELTA - inset
    floor_y = inset if inset_floor else 0.0

    points = []
    points.append([-half_floor, floor_y])
    points.append([half_floor, floor_y])
    for i in range(1, n_wall + 1):
        frac = i / n_wall
        y = floor_y + (wall_h - floor_y) * frac if inset_floor else frac * wall_h
        x = half_floor + (half_w - half_floor) * frac
        points.append([x, y])
    for i in range(1, n_arch):
        angle = np.pi * i / n_arch
        x = a * np.cos(angle)
        y = wall_h + b * np.sin(angle)
        points.append([x, y])
    for i in range(n_wall, 0, -1):
        frac = i / n_wall
        y = floor_y + (wall_h - floor_y) * frac if inset_floor else frac * wall_h
        x = half_floor + (half_w - half_floor) * frac
        points.append([-x, y])

    points = np.array(points)

    # Centre on cross-section centroid
    rect_area = TUNNEL_ALPHA * TUNNEL_WALL_HEIGHT
    rect_cy = TUNNEL_WALL_HEIGHT / 2
    a0 = TUNNEL_GAMMA / 2
    b0 = TUNNEL_DELTA
    ellipse_area = np.pi * a0 * b0 / 2
    ellipse_cy = TUNNEL_WALL_HEIGHT + 4 * b0 / (3 * np.pi)
    total_area = rect_area + ellipse_area
    centroid_y = (rect_area * rect_cy + ellipse_area * ellipse_cy) / total_area
    points[:, 1] -= centroid_y

    return points


# ============================================================
# Mesh extrusion along a 3D path
# ============================================================

def create_profile_mesh(path_points, profile_2d):
    """
    Extrude a 2D profile along a 3D path to create a triangulated mesh.

    Uses CMS convention: Y = up for the world-up reference vector.
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

        # CMS convention: Y is up
        if abs(tangent[1]) < 0.9:
            world_up = np.array([0, 1, 0])
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
# Tunnel centerline vertices (survey data, mm)
# ============================================================

_correctedVert = [
    (-86.57954338701529,    0.1882163986665546),
    (-1731.590867740335,    3.764327973349282),
    (-3549.761278867689,    7.716872345365118),
    (-5887.408950317142,    12.798715109387558),
    (-8053.403266181902,    -504.23173203003535),
    (-10046.991360867298,   -1282.5065405198511),
    (-11783.350377373874,   -2930.9057600491833),
    (-12913.652590171332,   -4580.622494369192),
    (-13095.344153684957,   -7536.749251839814),
    (-13099.610392054752,   -9015.000846973791),
    (-13278.792403586143,   -11101.567842600896),
    (-13372.39869252341,    -13536.146959364076),
    (-13292.093029091975,   -15710.234580371536),
    (-12779.140603923677,   -17972.21925955668),
    (-11659.12755425337,    -19887.69754879509),
    (-10105.714877251532,   -21630.204967658145),
    (-7512.845769209047,    -23201.0590309365),
    (-5262.530506741277,    -23466.820585854904),
    (-2751.72374851779,     -23472.278861416264),
    (-241.41890069074725,   -23651.64908934632),
    (1749.6596420124115,    -23742.93404270002),
    (3827.568683300815,     -23747.45123626804),
    (6078.6368113632525,    -23752.344862633392),
    (8502.613071001502,     -23844.570897980426),
    (11446.568501358292,    -23764.01427935077),
    (13438.399909656131,    -23594.431304151418),
    (15777.051401898476,    -23251.689242178036),
    (18289.614846509525,    -22648.455684448927),
    (20889.761655300477,    -21697.58643838109),
    (23143.841245741598,    -20659.00835053422),
    (25486.006110759066,    -19098.88262197991),
    (27742.09334278597,     -17364.656724658227),
    (28871.391734790544,    -16062.763895075637),
    (30781.662703665817,    -14153.873179790575),
    (32518.021720172394,    -12505.473960261239),
    (34513.49197884447,     -11075.029330388788),
    (36636.57295581305,     -10427.47081077351),
    (38759.40297758341,     -9866.868267342572),
    (41357.416667189485,    -9655.12481884172),
    (43694.93886103982,     -9703.684649697909),
    (46379.03018363646,     -9666.041369964427),
    (49409.43967978114,     -9629.150955825604),
    (51660.88424064092,     -9503.610617914434),
    (54258.0195870532,      -9596.213086058811),
    (57028.564975437745,    -9602.236010816167),
    (59539.87364405768,     -9433.782334008818),
    (62050.42944708294,     -9526.196585754526),
]

# Shift origin to CMS IP5 and convert mm → m
_X_SHIFT = 11908.8279764855    # mm
_Y_SHIFT = 13591.106147774964  # mm

correctedVertWithShift = [
    ((x - _X_SHIFT) / 1000, (y + _Y_SHIFT) / 1000)
    for x, y in _correctedVert
]

# Vertical offset of tunnel above CMS IP (m)
Y_POSITION = 22  # m


# ============================================================
# Build the default fiducial-volume mesh
# ============================================================

def build_fiducial_mesh(y_position=Y_POSITION,
                        detector_thickness=DETECTOR_THICKNESS):
    """
    Construct the fiducial-volume trimesh for the GRENDEL tunnel.

    Parameters
    ----------
    y_position : float
        Vertical (Y) offset of the tunnel centreline above the IP (m).
    detector_thickness : float
        Inset from the tunnel wall defining the fiducial volume (m).

    Returns
    -------
    mesh : trimesh.Trimesh
        Closed triangulated mesh of the fiducial volume.
    path_3d : ndarray, shape (N, 3)
        3D tunnel centreline coordinates.
    """
    # CMS convention: X-Z plane for tunnel path, Y = vertical
    path_3d = np.array(
        [[x, y_position, z] for x, z in correctedVertWithShift])

    profile = tunnel_profile_points(inset=detector_thickness,
                                    inset_floor=False)
    verts, faces = create_profile_mesh(path_3d, profile)
    mesh = trimesh.Trimesh(vertices=verts, faces=faces)
    # The Frenet-frame extrusion in create_profile_mesh produces inconsistent
    # face winding along curved sections; without fix_normals, mesh.volume
    # collapses to roughly 1/3 of the true volume (signed integral cancels).
    mesh.fix_normals()

    print(f"Fiducial volume built:")
    print(f"  Volume:    {mesh.volume:.1f} m³")
    print(f"  Y offset:  {y_position} m")
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
                world_up = np.array([0., 1., 0.])
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
            world_up = np.array([0., 1., 0.])
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
            world_up = np.array([0., 1., 0.])
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
        wu = np.array([0., 1., 0.]) if abs(seg_hat[1]) < 0.9 else np.array([0., 0., 1.])
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
        wu = np.array([0., 1., 0.]) if abs(seg_hat[1]) < 0.9 else np.array([0., 0., 1.])
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
