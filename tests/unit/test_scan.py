"""The scan driver against a synthetic model: a beam of long-lived particles
aimed at the detector, two-body decays from a template bundle, and a
ModelSpec that keeps every policy at its default."""
import json

import numpy as np
import pandas as pd
import pytest

from grendel.geometry.raycast import get_mesh
from grendel.reco.acceptance import build_event_mc
from grendel.reco.exclusion import find_exclusion_band
from grendel.scan import CouplingGrid, MassPoint, ModelSpec, ScanConfig, run_point, run_scan, secondary_threshold_columns

MASS = 1.0
FIELDS = ("u2_min", "u2_max", "u2_min_open", "u2_max_open", "peak_u2", "has_sensitivity")


def aimed_beam(path, n=3000, total_pb=3e7, seed=1):
    """A headerless weight,E,px,py,pz CSV of particles heading for the gallery."""
    rng = np.random.default_rng(seed)
    p = rng.uniform(5.0, 60.0, n)
    eta = rng.uniform(-2.5, 2.5, n)
    phi = np.pi / 2 + rng.uniform(-0.6, 0.6, n)
    pt = p / np.cosh(eta)
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    energy = np.sqrt(p**2 + MASS**2)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, np.column_stack([np.full(n, total_pb / n), energy, px, py, pz]), delimiter=",", fmt="%.8e")


def two_body_bundle(n_templates=200, ctau=1e-9, seed=2):
    """Isotropic rest-frame X -> mu pi decays in the common bundle format."""
    rng = np.random.default_rng(seed)
    m_mu, m_pi = 0.10566, 0.13957
    pstar = np.sqrt((MASS**2 - (m_mu + m_pi)**2) * (MASS**2 - (m_mu - m_pi)**2)) / (2 * MASS)
    cos = rng.uniform(-1, 1, n_templates)
    azimuth = rng.uniform(0, 2 * np.pi, n_templates)
    sin = np.sqrt(1 - cos**2)
    d = np.column_stack([sin * np.cos(azimuth), sin * np.sin(azimuth), cos]) * pstar
    mom = np.empty((2 * n_templates, 3))
    mom[0::2], mom[1::2] = d, -d
    e = np.empty(2 * n_templates)
    e[0::2], e[1::2] = np.sqrt(pstar**2 + m_mu**2), np.sqrt(pstar**2 + m_pi**2)
    return {"daughter_counts": np.full(n_templates, 2, np.int32), "pdg": np.tile([13, 211], n_templates).astype(np.int32),
            "px": mom[:, 0], "py": mom[:, 1], "pz": mom[:, 2], "energy": e,
            "mass": np.tile([m_mu, m_pi], n_templates), "charge": np.tile([-1.0, 1.0], n_templates),
            "stable": np.ones(2 * n_templates, bool), "mass_GeV": np.float64(MASS), "ctau_m_u2eq1": np.float64(ctau)}


class TemplateDecays:
    name = "synthetic-two-body"

    def __init__(self, bundle):
        self.bundle = bundle
        self.ctau_ref = float(bundle["ctau_m_u2eq1"])

    def build(self, p4, direction, entry_d, exit_d, n_samples, rng, return_mc=False):
        out = build_event_mc(p4, direction, entry_d, exit_d, self.bundle, n_samples, rng, return_mc=return_mc)
        return out if return_mc else (*out, None)

    def sample_weights(self, template_index):
        return None


class ToySpec(ModelSpec):
    name = "toy"
    grid = CouplingGrid(-11.0, -3.0, 80)

    def __init__(self, root, bundle, chunk=None, seed=7):
        self.root, self.bundle, self.chunk, self.seed = root, bundle, chunk, seed

    def vectors_path(self, pt):
        return self.root / "vectors" / f"m_{pt.label}.csv"

    def geometry_cache_path(self, pt, vectors_path):
        return self.root / "geometry" / f"geom_{pt.label}.npz"

    def decay_backend(self, pt, vectors_path):
        return TemplateDecays(self.bundle)

    def rng_for(self, pt, cfg):
        return np.random.default_rng([self.seed, hash(cfg.seed_salt) % 2**31])

    def chunks(self, pt, n_events, cfg, rng):
        if not self.chunk:
            return [(slice(0, n_events), rng)]
        return [(slice(s, min(s + self.chunk, n_events)), rng) for s in range(0, n_events, self.chunk)]

    def base_row(self, pt, n_events, n_hits):
        return {"mass_GeV": pt.mass, "n_events": n_events, "n_hits": n_hits}

    def empty_row(self, pt, n_events, n_hits, cfg, *, evaluated=False):
        return {"mass_GeV": pt.mass, "n_events": n_events, "n_hits": n_hits, "has_sensitivity": False,
                "peak_N": float("nan")}

    def finish(self, pt, base, backend, arrays, cfg):
        result = find_exclusion_band(arrays.grid, arrays.N, cfg.primary)
        secondary_threshold_columns(result, lambda t: find_exclusion_band(arrays.grid, arrays.N, t),
                                    cfg.secondary, FIELDS)
        result.update(base)
        return result


@pytest.fixture(scope="module")
def mesh():
    return get_mesh()


@pytest.fixture
def toy(tmp_path):
    spec = ToySpec(tmp_path, two_body_bundle())
    aimed_beam(spec.vectors_path(MassPoint(MASS)))
    return spec


def test_run_point_solves_a_nested_pair_of_islands(toy, mesh):
    cfg = ScanConfig(decay_samples=5, thresholds=(3.0, 10.0))
    result = run_point(toy, MassPoint(MASS), cfg, mesh, keep_arrays=True)
    row = result.row
    assert row["has_sensitivity"] and row["n_hits"] > 0 and row["n_events"] == 3000
    assert row["peak_N"] > 10 and 0 < row["u2_min"] < row["peak_u2"] < row["u2_max"]
    assert row["u2_min"] <= row["u2_min_N10"] < row["u2_max_N10"] <= row["u2_max"]
    assert toy.geometry_cache_path(MassPoint(MASS), None).is_file()
    arrays = result.arrays
    assert arrays.N.shape == (80,) and arrays.d.shape == (row["n_hits"], 5)
    assert arrays.hit_estimator == "exact" and arrays.n_hits_eval == row["n_hits"]


def test_the_same_seed_reproduces_and_a_salt_changes_the_draw(toy, mesh):
    cfg = ScanConfig(decay_samples=5)
    first = run_point(toy, MassPoint(MASS), cfg, mesh).row
    again = run_point(toy, MassPoint(MASS), cfg, mesh).row
    salted = run_point(toy, MassPoint(MASS), ScanConfig(decay_samples=5, seed_salt="control-1"), mesh).row
    assert first == again
    assert salted["peak_N"] != first["peak_N"]
    assert salted["u2_min"] == pytest.approx(first["u2_min"], rel=0.3)   # same physics, other statistics


def test_accumulate_mode_sums_the_chunk_yields_exactly(tmp_path, mesh):
    # A model that chunks its events draws different random numbers from a
    # one-pass model (the acceptance Monte Carlo consumes the generator per
    # batch), which is why chunking is part of each model's numerical
    # policy. What the driver guarantees is that the per-chunk yield curves
    # it accumulates equal one scan over the concatenated arrays.
    from grendel.constants import L_INT_PB
    from grendel.reco.acceptance import scan_u2
    chunked = ToySpec(tmp_path, two_body_bundle(), chunk=17)
    chunked.scan_mode = "accumulate"
    aimed_beam(chunked.vectors_path(MassPoint(MASS)))
    arrays = run_point(chunked, MassPoint(MASS), ScanConfig(decay_samples=4), mesh, keep_arrays=True).arrays
    _, one_pass = scan_u2(arrays.d, arrays.passed, arrays.path, arrays.weights, arrays.beta_gamma,
                          arrays.ctau_ref, L_INT_PB, arrays.grid)
    assert np.allclose(arrays.N, one_pass, rtol=1e-13, atol=0.0)
    assert arrays.d.shape[0] == arrays.n_hits_eval > 17


def test_missing_vectors_skip_the_point_and_are_reported(toy, mesh, tmp_path):
    assert run_point(toy, MassPoint(2.0), ScanConfig(decay_samples=3), mesh) is None
    out = run_scan(toy, [MassPoint(MASS), MassPoint(2.0)], ScanConfig(decay_samples=3), tmp_path / "out",
                   verbose=False)
    meta = json.loads((out.parent / "run_metadata.json").read_text())
    assert meta["n_points_processed"] == 1 and meta["n_points_skipped"] == 1
    assert meta["skipped"][0]["reason"] == "missing_or_empty_vectors"


def test_run_scan_writes_the_layout_and_resumes(toy, mesh, tmp_path):
    cfg = ScanConfig(decay_samples=3, thresholds=(3.0, 10.0))
    out = run_scan(toy, [MassPoint(MASS)], cfg, tmp_path / "out", verbose=False)
    assert out.name == "sensitivity.csv"
    assert {p.name for p in out.parent.iterdir()} >= {"sensitivity.csv", "run_metadata.json", "scan_status.json"}
    frame = pd.read_csv(out)
    assert list(frame["mass_GeV"]) == [MASS] and bool(frame.loc[0, "has_sensitivity"])
    assert {"u2_min_N10", "has_sensitivity_N10"} <= set(frame.columns)
    meta = json.loads((out.parent / "run_metadata.json").read_text())
    assert meta["model"] == "toy" and meta["thresholds"] == [3.0, 10.0] and meta["n_sensitive"] == 1
    # Resuming keeps the row instead of recomputing it.
    (tmp_path / "out" / "geometry_cache").mkdir(exist_ok=True)
    toy.geometry_cache_path(MassPoint(MASS), None).unlink()   # a rerun would have to ray-cast again
    run_scan(toy, [MassPoint(MASS)], cfg, tmp_path / "out", resume=True, verbose=False)
    resumed = json.loads((out.parent / "run_metadata.json").read_text())
    assert resumed["n_points_requested"] == resumed["n_points_processed"] == 1
    assert not toy.geometry_cache_path(MassPoint(MASS), None).exists()
    assert pd.read_csv(out).equals(frame)
