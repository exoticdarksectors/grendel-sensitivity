import numpy as np
import pytest

from grendel.constants import CMS_ORIGIN
from grendel.geometry import raycast
from grendel.io.vectors import format_mass_for_filename, parse_mass_from_filename


def test_mesh_is_a_closed_volume():
    mesh = raycast.get_mesh()
    assert mesh.is_watertight
    assert 600 < mesh.volume < 800          # m^3, the fiducial air volume


def test_directions_are_unit_vectors():
    rng = np.random.default_rng(1)
    eta = rng.uniform(-4, 4, 1000)
    phi = rng.uniform(-np.pi, np.pi, 1000)
    d = raycast.directions_from_eta_phi(eta, phi)
    assert d.shape == (1000, 3)
    assert np.allclose(np.linalg.norm(d, axis=1), 1.0)
    # eta = 0 is transverse, large eta is along +z
    assert abs(raycast.directions_from_eta_phi(np.zeros(1), np.zeros(1))[0, 2]) < 1e-12
    assert raycast.directions_from_eta_phi(np.array([6.0]), np.zeros(1))[0, 2] > 0.99


def test_raycast_finds_the_detector_above_the_ip():
    mesh = raycast.get_mesh()
    # the gallery sits 23.5-25 m above IP5 but off the vertical: rays going
    # straight up miss it, rays going down cannot hit it; eta = 0.44 at
    # phi = 1.50 points into it
    up = raycast.compute_geometry(np.array([0.0, 0.44]), np.array([np.pi / 2, 1.50]), mesh, CMS_ORIGIN)
    down = raycast.compute_geometry(np.zeros(3), np.full(3, -np.pi / 2), mesh, CMS_ORIGIN)
    assert not down[0].any()
    assert list(up[0]) == [False, True]
    assert 10.0 < up[1][1, 0] < up[2][1, 0]


def _two_boxes():
    import trimesh
    return trimesh.util.concatenate([trimesh.creation.box((2.0, 2.0, 2.0)).apply_translation([0.0, 10.0, 0.0]),
                                     trimesh.creation.box((2.0, 3.0, 2.0)).apply_translation([0.0, 15.0, 0.0])])


def test_a_ray_crossing_the_volume_twice_keeps_both_segments():
    # the as-built gallery is crossed twice by almost no ray from IP5 (one in
    # 40 000 directions, by under 0.5 m), so the logic is tested on two boxes
    # straight above the IP
    hits, entry_d, exit_d = raycast.compute_geometry(np.zeros(1), np.full(1, np.pi / 2), _two_boxes(), CMS_ORIGIN)
    assert hits[0] and entry_d.shape == (1, 2)
    assert np.allclose([entry_d[0, 0], exit_d[0, 0], entry_d[0, 1], exit_d[0, 1]],
                       [9.0, 11.0, 13.5, 16.5], atol=1e-9)
    assert raycast.path_length(entry_d, exit_d)[0] == pytest.approx(5.0, abs=1e-9)


def test_decay_distances_cover_every_segment_and_keep_single_segment_draws():
    entry = np.array([[10.0, np.nan], [10.0, 13.0]])
    exit_ = np.array([[12.0, np.nan], [12.0, 19.0]])
    d = raycast.sample_decay_distances(entry, exit_, 20000, np.random.default_rng(5))
    legacy = np.random.default_rng(5).uniform(entry[:, :1], exit_[:, :1], size=(2, 20000))
    assert np.array_equal(d[0], legacy[0])                     # one segment: the historical draw
    two = d[1]
    assert np.all(((two >= 10) & (two < 12)) | ((two >= 13) & (two < 19)))
    assert np.mean(two < 12) == pytest.approx(2 / 8, abs=0.01)  # uniform over 8 m of path
    assert np.allclose(raycast.path_length(entry, exit_), [2.0, 8.0])


def test_the_reconstruction_tests_containment_in_the_ray_cast_mesh():
    from grendel.geometry import fiducial, reco_common
    mesh = raycast.get_mesh()
    assert reco_common.points_in_fiducial is fiducial.points_in_fiducial
    rng = np.random.default_rng(4)
    pts = rng.uniform(*mesh.bounds, size=(4000, 3))
    assert np.array_equal(fiducial.points_in_fiducial(pts), mesh.contains(pts))
    # 19.75 cm outside the mesh at a bend, accepted by the analytic test
    outside = np.array([[-26.199238, 20.882130, 9.127790]])
    assert not fiducial.points_in_fiducial(outside)[0] and not mesh.contains(outside)[0]


def test_containment_counts_a_crossing_through_a_shared_edge_or_vertex_once():
    from grendel.geometry.fiducial import VerticalRayParity
    mesh = raycast.get_mesh()
    test = VerticalRayParity(mesh)
    v = np.asarray(mesh.vertices)
    edges = np.sort(mesh.edges_unique, axis=1)
    # 1 m below every vertex and every edge midpoint: the vertical ray from
    # the point runs through that vertex or edge, shared by several triangles
    pts = np.vstack([v, 0.5 * (v[edges[:, 0]] + v[edges[:, 1]])]) - [0.0, 1.0, 0.0]
    nudges = np.array([[1e-6, 0, 1.3e-6], [-1e-6, 0, 1.3e-6], [1e-6, 0, -1.3e-6], [-1e-6, 0, -1.3e-6]])
    ref = np.array([mesh.contains(pts + d) for d in nudges])
    clear = ref.min(axis=0) == ref.max(axis=0)       # off the surface: every nudge agrees
    assert clear.sum() > 0.9 * len(pts) and ref[0, clear].sum() > 1000
    assert np.array_equal(test(pts[clear]), ref[0, clear])


def test_geometry_cache_roundtrip_is_exact(tmp_path):
    mesh = raycast.get_mesh()
    rng = np.random.default_rng(2)
    eta = rng.uniform(-2, 2, 2000)
    phi = rng.uniform(0, np.pi, 2000)
    cache = tmp_path / "geom.npz"
    first = raycast.load_or_compute_geometry(cache, eta, phi, mesh)
    assert cache.exists()
    again = raycast.load_or_compute_geometry(cache, eta, phi, mesh)
    for a, b in zip(first, again):
        assert np.array_equal(a, b, equal_nan=True)
    # a newer source invalidates the cache
    stale = raycast.load_or_compute_geometry(cache, eta, phi, mesh,
                                             source_mtime=cache.stat().st_mtime + 10)
    for a, b in zip(first, stale):
        assert np.array_equal(a, b, equal_nan=True)
    # other directions of the same length are not served the cached result
    other = raycast.load_or_compute_geometry(cache, eta[::-1].copy(), phi[::-1].copy(), mesh)
    assert np.array_equal(other[0], first[0][::-1])


def test_a_cache_that_cannot_be_verified_is_recomputed(tmp_path):
    mesh = _two_boxes()
    eta = np.array([0.0, 0.0])
    phi = np.array([-1.0, np.pi / 2])
    truth = raycast.compute_geometry(eta, phi, mesh, CMS_ORIGIN)
    assert list(truth[0]) == [False, True]
    # A cache from before the fingerprint, holding misses for both rays (as
    # one cast for other directions would): it is not trusted.
    cache = tmp_path / "old.npz"
    np.savez_compressed(cache, hits=np.zeros(2, bool), entry_d=np.full(2, np.nan), exit_d=np.full(2, np.nan))
    loaded = raycast.load_or_compute_geometry(cache, eta, phi, mesh)
    for a, b in zip(loaded, truth):
        assert np.array_equal(a, b, equal_nan=True)
    assert "fingerprint" in np.load(cache).files and np.load(cache)["entry_d"].shape == (2, 2)


def test_a_cache_for_another_origin_or_mesh_is_not_served(tmp_path):
    mesh = raycast.get_mesh()
    eta, phi = np.array([0.44]), np.array([1.50])
    cache = tmp_path / "geom.npz"
    first = raycast.load_or_compute_geometry(cache, eta, phi, mesh)
    moved = np.asarray(CMS_ORIGIN, dtype=float) + np.array([0.0, 1.0, 0.0])
    shifted = raycast.load_or_compute_geometry(cache, eta, phi, mesh, origin=moved)
    assert shifted[1][0, 0] != first[1][0, 0]
    assert np.array_equal(shifted[1], raycast.compute_geometry(eta, phi, mesh, moved)[1])
    edited = mesh.copy()
    edited.apply_translation([0.0, 0.5, 0.0])
    other = raycast.load_or_compute_geometry(cache, eta, phi, edited, origin=moved)
    assert np.array_equal(other[1], raycast.compute_geometry(eta, phi, edited, moved)[1])
    assert other[1][0, 0] != shifted[1][0, 0]


@pytest.mark.parametrize("mass,label", [(0.2, "0p200"), (0.305, "0p305"),
                                        (1.025, "1p025"), (3.62, "3p620"), (10.0, "10p000")])
def test_mass_labels_roundtrip(mass, label):
    assert format_mass_for_filename(mass) == label
    assert parse_mass_from_filename(label) == mass
    assert parse_mass_from_filename(f"mN_{label}") == mass
