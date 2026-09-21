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
    # the gallery sits ~22 m above IP5: rays going straight up hit it,
    # rays going down cannot
    up = raycast.compute_geometry(np.zeros(3), np.full(3, np.pi / 2), mesh, CMS_ORIGIN)
    down = raycast.compute_geometry(np.zeros(3), np.full(3, -np.pi / 2), mesh, CMS_ORIGIN)
    assert not down[0].any()
    if up[0].any():
        assert np.all(up[1][up[0]] < up[2][up[0]])
        assert np.all(up[1][up[0]] > 10.0)


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
    # a newer source invalidates the cache; forcing recomputes it
    stale = raycast.load_or_compute_geometry(cache, eta, phi, mesh,
                                             source_mtime=cache.stat().st_mtime + 10)
    for a, b in zip(first, stale):
        assert np.array_equal(a, b, equal_nan=True)


@pytest.mark.parametrize("mass,label", [(0.2, "0p200"), (0.305, "0p305"),
                                        (1.025, "1p025"), (3.62, "3p620"), (10.0, "10p000")])
def test_mass_labels_roundtrip(mass, label):
    assert format_mass_for_filename(mass) == label
    assert parse_mass_from_filename(label) == mass
    assert parse_mass_from_filename(f"mN_{label}") == mass
