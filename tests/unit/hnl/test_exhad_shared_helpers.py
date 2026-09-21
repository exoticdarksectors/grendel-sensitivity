"""Shared exHad template helpers (no exHad runtime needed)."""
import numpy as np
import pytest

from grendel.reco import templates as xt


class _FlakyGenerator:
    """Fails the second sub-request once, like exHad's rare per-event consistency error."""

    def __init__(self):
        self.calls = []

    def generate_all(self, mass, events, *, seed, **kwargs):
        self.calls.append((mass, events, seed))
        if len(self.calls) == 2:
            raise RuntimeError("RuntimeError: Full decay does not conserve four-momentum")
        return {"events": [[[0.0, 0.0, 0.0, mass, mass, 22]] for _ in range(events)],
                "channel_labels": ["2gamma"] * events}


def test_generate_all_robust_retries_only_the_failed_sub_request():
    g = _FlakyGenerator()
    events, labels, retries = xt.generate_all_robust(g, 1.5, 5000, seed=7, sub_request=2000, log=lambda *_: None)
    assert len(events) == 5000 and len(labels) == 5000
    assert [c[1] for c in g.calls] == [2000, 2000, 2000, 1000]
    assert len(retries) == 1 and retries[0][:2] == (1, 0)
    # sub-request seeds are distinct and change on retry
    seeds = [c[2] for c in g.calls]
    assert len(set(seeds)) == 4
    assert g.calls[2][2] == xt.sub_request_seed(7, 1, 1)


def test_generate_all_robust_gives_up_after_max_retries():
    class Broken:
        def generate_all(self, *a, **k):
            raise RuntimeError("always")
    with pytest.raises(RuntimeError, match="failed 2 times"):
        xt.generate_all_robust(Broken(), 1.0, 10, seed=1, max_retries=1, log=lambda *_: None)


def test_point_seed_namespaces_are_independent_and_stable():
    a = xt.point_seed(1234, "Ue", "1p000", namespace="exhad-hnl-templates/v1")
    b = xt.point_seed(1234, "Ue", "1p000", namespace="exhad-bc4-templates/v1")
    assert a != b and 0 <= a < 2 ** 63
    assert a == xt.point_seed(1234, "Ue", "1p000", namespace="exhad-hnl-templates/v1")


def test_loglog_interp_is_exact_on_power_laws_and_refuses_extrapolation():
    xs = np.array([1.0, 4.0]); ys = np.array([2.0, 32.0])  # y = 2 x^2
    assert xt.loglog_interp(2.0, xs, ys) == pytest.approx(2.0 * 2.0 ** 2)
    with pytest.raises(ValueError):
        xt.loglog_interp(8.0, xs, ys)


def test_channel_fractions_sorted_descending():
    assert list(xt.channel_fractions(["a", "b", "b", "c", "b"]).items())[0] == ("b", 0.6)


def _bundle(path, mass, seed, labels, charges, flavor="BC4"):
    n = len(charges)
    np.savez_compressed(
        path, daughter_counts=np.array([n], np.int32), pdg=np.array([211] * n, np.int32),
        px=np.zeros(n), py=np.zeros(n), pz=np.zeros(n), energy=np.ones(n), mass=np.zeros(n),
        charge=np.array(charges, float), stable=np.ones(n, bool),
        channel_label=np.array(labels, dtype=str), n_charged=np.array([int(sum(abs(q) > 0.5 for q in charges))], np.int16),
        mass_GeV=np.float64(mass), ctau_m_u2eq1=np.float64(1e-3), ctau_m_u2eq1_reference=np.float64(2e-3),
        exhad_seed=np.uint64(seed), flavor=np.array(flavor))


def test_complete_manifest_points_merges_prior_records_and_rebuilds_the_rest(tmp_path):
    _bundle(tmp_path / "templates_0p500.npz", 0.5, 11, ["pipi"], [1.0, -1.0])
    _bundle(tmp_path / "templates_1p000.npz", 1.0, 12, ["KK"], [0.0, 0.0])
    _bundle(tmp_path / "templates_2p000.npz", 2.0, 13, ["hadronic"], [1.0, -1.0, 0.0])
    prior = [{"mass_GeV": 0.5, "file": "templates_0p500.npz", "exhad_seed": 11, "retries": [[0, 0, "x"]]},
             {"mass_GeV": 9.0, "file": "templates_9p000.npz", "error": "failed"}]
    generated = [{"mass_GeV": 2.0, "file": "templates_2p000.npz", "exhad_seed": 13, "retries": []}]
    points, counts = xt.complete_manifest_points(tmp_path, generated, prior)
    assert counts == {"generated": 1, "from_prior_manifest": 1, "rebuilt_from_bundle": 1}
    assert [p["file"] for p in points] == ["templates_0p500.npz", "templates_1p000.npz", "templates_2p000.npz"]
    assert points[0]["retries"] == [[0, 0, "x"]]          # carried over, with its retry log
    rebuilt = points[1]
    assert rebuilt["exhad_seed"] == 12 and rebuilt["vis_frac_ge2_charged"] == 0.0
    assert rebuilt["ctau_m_u2eq1_reference"] == 2e-3 and rebuilt["channel_fractions"] == {"KK": 1.0}
    assert "rebuilt" in rebuilt["record_source"]
    assert "flavor" not in rebuilt                        # BC4 bundles carry no HNL flavour


def test_record_from_bundle_keeps_the_hnl_flavour_and_relative_path(tmp_path):
    (tmp_path / "Umu").mkdir()
    _bundle(tmp_path / "Umu" / "templates_1p000.npz", 1.0, 5, ["Pimu"], [1.0, -1.0], flavor="Umu")
    rec = xt.record_from_bundle(tmp_path / "Umu" / "templates_1p000.npz", tmp_path)
    assert rec["flavor"] == "Umu" and rec["file"] == "Umu/templates_1p000.npz"
    assert rec["vis_frac_ge2_charged"] == 1.0 and rec["mean_multiplicity"] == 2.0


def test_generation_runs_reconstructs_a_legacy_manifest_entry():
    legacy = {"generated_at": "t0", "points": [1, 2, 3], "max_retries": 10, "elapsed_s": 5.0,
              "seed_policy": "...; sub-requests of 250 events seeded by ..."}
    runs = xt.generation_runs(legacy, {"generated_at": "t1", "n_generated": 1})
    assert runs[0]["n_generated"] == 3 and runs[0]["sub_request"] == 250 and runs[0]["max_retries"] == 10
    assert runs[1] == {"generated_at": "t1", "n_generated": 1}
    assert xt.generation_runs({}, {"generated_at": "t1"}) == [{"generated_at": "t1"}]
    kept = {"generation_runs": [{"generated_at": "a"}], "points": []}
    assert xt.generation_runs(kept, {"generated_at": "b"}) == [{"generated_at": "a"}, {"generated_at": "b"}]
