"""The BC4 variation campaign: variation discovery, seeds, the single-source
envelope rule, and the retained-run bookkeeping."""
import hashlib
import json
import os
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from grendel.band.campaign_store import git_head, stable_seed
from grendel.band.scalar import campaign as uncertainty_band
from grendel.band.scalar.campaign import (CODE_INPUTS, CONFIG_KEYS, campaign_config, config_keys,
                                          finalize_and_compact, load_recorded_config, templates_provenance,
                                          validate_retained_completion)
from grendel.band.scalar.combine import combine_band
from grendel.band.scalar.variations import (DECAY_VARIATION, EXHAD_DECAY_SCHEME, EXHAD_DECAY_VARIATION,
                                            NUMERICAL_CONTROL_VARIATIONS, discover_variations)
from grendel.io.paths import repo_root


@pytest.mark.external
def test_complete_fonll_manifest_is_discovered():
    grid_dir = os.environ.get("GRENDEL_FONLL_GRID_DIR")
    if not grid_dir or not (Path(grid_dir) / "variation_manifest.json").is_file():
        pytest.skip("GRENDEL_FONLL_GRID_DIR does not hold a complete FONLL variation campaign")
    variations = discover_variations(Path(grid_dir))
    counts = pd.Series([variation["axis"] for variation in variations]).value_counts()
    assert counts.to_dict() == {"pdf": 100, "scale": 6, "mass": 2, "numerical_control": 2, "central": 1,
                                "decay_model": 1}
    assert variations[0]["name"] == "central"
    assert variations[1]["name"] == DECAY_VARIATION
    assert variations[0]["path"] == variations[1]["path"]
    assert variations[0]["width_scheme"] != variations[1]["width_scheme"]
    assert tuple(v["name"] for v in variations[-2:]) == NUMERICAL_CONTROL_VARIATIONS


def test_variation_seeds_are_deterministic_and_independent():
    central = stable_seed(42, "central", "parent_pool")
    assert central == stable_seed(42, "central", "parent_pool")
    assert central != stable_seed(42, DECAY_VARIATION, "parent_pool")
    assert central != stable_seed(42, NUMERICAL_CONTROL_VARIATIONS[0], "parent_pool")
    assert stable_seed(42, NUMERICAL_CONTROL_VARIATIONS[0], "parent_pool") != \
        stable_seed(42, NUMERICAL_CONTROL_VARIATIONS[1], "parent_pool")
    assert central != stable_seed(42, "central", "0p500", "production")
    assert 0 <= central < 2**32


def _row(name, axis, umin, umax):
    return {"mass_GeV": 1.0, "variation": name, "axis": axis, "has_sensitivity": True,
            "u2_min": umin, "u2_max": umax, "u2_min_open": False, "u2_max_open": False}


def _reference(mass=1.0, sensitive=True, u2_min=2e-8, u2_max=2e-4):
    return pd.DataFrame([{"mass_GeV": mass, "has_sensitivity": sensitive,
                          "u2_min": u2_min if sensitive else float("nan"),
                          "u2_max": u2_max if sensitive else float("nan"),
                          "u2_min_open": False, "u2_max_open": False}])


def test_combine_band_rebases_single_source_intervals():
    raw = pd.DataFrame([
        _row("central", "central", 1e-8, 1e-4),
        _row("scale_up", "scale", 10**-7.8, 10**-3.8),
        _row("scale_down", "scale", 10**-8.1, 10**-4.1),
        _row("pdf_1", "pdf", 10**-7.95, 10**-3.95),
        _row("pdf_2", "pdf", 10**-8.05, 10**-4.05),
        _row("mb_dn", "mass", 10**-7.9, 10**-3.9),
        _row("mb_up", "mass", 10**-8.08, 10**-4.08),
        _row(DECAY_VARIATION, "decay_model", 10**-8.3, 10**-3.7),
    ])
    out = combine_band(raw, _reference()).iloc[0]
    assert out["envelope_definition"] == "single_source_variation_envelope"
    assert bool(out["any_variation_sensitive"])
    assert out["u2_min_central"] == pytest.approx(2e-8)
    assert out["u2_min_campaign_central"] == pytest.approx(1e-8)
    assert out["u2_min_scale_envelope_lo_dex"] == pytest.approx(-0.1)
    assert out["u2_min_scale_envelope_hi_dex"] == pytest.approx(0.2)
    assert out["u2_min_pdf_p16_shift_dex"] == pytest.approx(-0.034)
    assert out["u2_min_pdf_p84_shift_dex"] == pytest.approx(0.034)
    assert out["u2_min_pdf_std_dex"] == pytest.approx(0.05 * 2**0.5)
    assert out["u2_min_bottom_mass_envelope_lo_dex"] == pytest.approx(-0.08)
    assert out["u2_min_bottom_mass_envelope_hi_dex"] == pytest.approx(0.1)
    assert out["u2_min_decay_model_shift_dex"] == pytest.approx(-0.3)
    assert out["u2_min_decay_model_alt"] == pytest.approx(2e-8 * 10**-0.3)
    assert out["u2_min_envelope_lo"] == pytest.approx(2e-8 * 10**-0.3)
    assert out["u2_min_envelope_hi"] == pytest.approx(2e-8 * 10**0.2)
    assert out["u2_min_envelope_lo_source"] == f"decay_model:{DECAY_VARIATION}"
    assert out["u2_min_envelope_hi_source"] == "scale:scale_up"
    assert out["u2_max_envelope_hi"] == pytest.approx(2e-4 * 10**0.3)
    assert out["u2_max_envelope_hi_source"] == f"decay_model:{DECAY_VARIATION}"


def test_pdf_replica_outlier_does_not_define_display_envelope():
    rows = [_row("central", "central", 1e-8, 1e-4)]
    rows.extend(_row(f"pdf_{index}", "pdf", 1e-8, 1e-4) for index in range(99))
    rows.append(_row("pdf_99", "pdf", 1e-5, 1e-1))
    rows.append(_row(DECAY_VARIATION, "decay_model", 1e-8, 1e-4))
    rows.append(_row(NUMERICAL_CONTROL_VARIATIONS[0], "numerical_control", 1e-5, 1e-1))
    rows.append(_row(NUMERICAL_CONTROL_VARIATIONS[1], "numerical_control", 1e-11, 1e-7))
    out = combine_band(pd.DataFrame(rows), _reference()).iloc[0]
    assert out["u2_min_pdf_p84"] == pytest.approx(2e-8)
    assert out["u2_min_envelope_hi"] == pytest.approx(2e-8)
    assert out["u2_min_pdf_std_dex"] > 0
    assert out["u2_min_numerical_repeat_median_abs_dex"] == pytest.approx(3.0)
    assert out["u2_min_numerical_repeat_max_abs_dex"] == pytest.approx(3.0)
    assert bool(out["u2_min_numerical_repeat_not_subdominant"])


def test_numerical_repeat_does_not_restore_physics_topology():
    raw = pd.DataFrame([
        {**_row("central", "central", 1e-8, 1e-4), "has_sensitivity": False},
        {**_row(DECAY_VARIATION, "decay_model", 1e-8, 1e-4), "has_sensitivity": False},
        _row(NUMERICAL_CONTROL_VARIATIONS[0], "numerical_control", 1e-8, 1e-4),
        {**_row(NUMERICAL_CONTROL_VARIATIONS[1], "numerical_control", 1e-8, 1e-4), "has_sensitivity": False},
    ])
    out = combine_band(raw, _reference(sensitive=False)).iloc[0]
    assert not bool(out["any_variation_sensitive"])
    assert bool(out["numerical_control_any_sensitive"])
    assert not bool(out["numerical_control_all_sensitive"])
    assert out["numerical_control_n_sensitive"] == 1
    assert out["numerical_control_sensitive_variations"] == NUMERICAL_CONTROL_VARIATIONS[0]
    assert bool(out["numerical_control_topology_differs_from_campaign"])
    assert bool(out["numerical_control_topology_differs_from_canonical"])
    assert out["numerical_control_topology_difference_variations_from_campaign"] == NUMERICAL_CONTROL_VARIATIONS[0]
    assert out["numerical_control_topology_difference_variations_from_canonical"] == NUMERICAL_CONTROL_VARIATIONS[0]
    assert out[f"u2_min_{NUMERICAL_CONTROL_VARIATIONS[0]}"] == pytest.approx(1e-8)
    assert not bool(out["numerical_control_included_in_envelope"])


def test_physical_variation_topology_change_is_named_at_canonical_gap():
    raw = pd.DataFrame([{**_row("central", "central", 1e-8, 1e-4), "has_sensitivity": False},
                        _row(DECAY_VARIATION, "decay_model", 1e-8, 1e-4)])
    raw["mass_GeV"] = 3.8
    out = combine_band(raw, _reference(mass=3.8, sensitive=False)).iloc[0]
    assert bool(out["any_variation_sensitive"])
    assert bool(out["physical_variation_any_sensitive"])
    assert out["physical_variation_n_sensitive"] == 1
    assert out["physical_variation_sensitive_variations"] == DECAY_VARIATION
    assert bool(out["physical_variation_topology_differs_from_campaign"])
    assert bool(out["physical_variation_topology_differs_from_canonical"])
    assert out["physical_variation_topology_difference_variations_from_campaign"] == DECAY_VARIATION
    assert out["physical_variation_topology_difference_variations_from_canonical"] == DECAY_VARIATION
    assert pd.isna(out["u2_min_envelope_lo"])


def test_combine_band_spans_every_decay_model_alternate():
    # Two decay-model alternates: the interval is the span of central and both,
    # the plain columns still describe the first alternate, and each alternate
    # gets its own named columns.
    raw = pd.DataFrame([
        _row("central", "central", 1e-8, 1e-4),
        _row(DECAY_VARIATION, "decay_model", 10**-8.3, 10**-3.7),
        _row(EXHAD_DECAY_VARIATION, "decay_model", 10**-7.9, 10**-4.4),
    ])
    out = combine_band(raw, _reference()).iloc[0]
    assert out["u2_min_decay_model_shift_dex"] == pytest.approx(-0.3)
    assert out[f"u2_min_decay_model_{DECAY_VARIATION}_shift_dex"] == pytest.approx(-0.3)
    assert out[f"u2_min_decay_model_{EXHAD_DECAY_VARIATION}_shift_dex"] == pytest.approx(0.1)
    assert out[f"u2_min_decay_model_{EXHAD_DECAY_VARIATION}"] == pytest.approx(2e-8 * 10**0.1)
    assert out["u2_min_decay_model_envelope_lo"] == pytest.approx(2e-8 * 10**-0.3)
    assert out["u2_min_decay_model_envelope_hi"] == pytest.approx(2e-8 * 10**0.1)
    assert out["u2_min_decay_model_envelope_lo_source"] == DECAY_VARIATION
    assert out["u2_min_decay_model_envelope_hi_source"] == EXHAD_DECAY_VARIATION
    assert out["u2_max_decay_model_envelope_lo"] == pytest.approx(2e-4 * 10**-0.4)
    assert out["u2_max_decay_model_envelope_lo_source"] == EXHAD_DECAY_VARIATION
    assert out["u2_max_decay_model_envelope_hi_source"] == DECAY_VARIATION
    assert out["u2_min_envelope_lo_source"] == f"decay_model:{DECAY_VARIATION}"
    assert out["u2_min_envelope_hi_source"] == f"decay_model:{EXHAD_DECAY_VARIATION}"
    assert out["u2_min_envelope_hi"] == pytest.approx(2e-8 * 10**0.1)
    assert out["physical_variation_sensitive_variations"] == f"{DECAY_VARIATION};{EXHAD_DECAY_VARIATION}"


# ------------------------------------------------------- retained runs ------

def test_completed_variation_is_checksummed_then_compacted(tmp_path):
    run_dir = tmp_path / "runs" / "central"
    vector_dir = run_dir / "llp_4vectors"
    geometry_dir = run_dir / "geometry_cache"
    result_dir = run_dir / "results"
    vector_dir.mkdir(parents=True)
    geometry_dir.mkdir()
    result_dir.mkdir()
    mass = 0.5
    config = {"config_sha256": "config-sha", "base_seed": 42, "parent_pool_seed": 123,
              "code_sha256": {"code.py": "abc"}, "grid_sha256": "grid-sha", "ray_backend": {"implementation": "test"}}
    vector = vector_dir / "mS_0p500.csv"
    vector.write_text("1,2,3,4,5\n")
    (vector_dir / "mS_0p500.meta.json").write_text(json.dumps({
        "config_sha256": config["config_sha256"], "mass_GeV": mass, "bytes": vector.stat().st_size,
        "sha256": hashlib.sha256(vector.read_bytes()).hexdigest()}))
    (geometry_dir / "geom_mS_0p500.npz").write_bytes(b"geometry")
    (result_dir / "mS_0p500.json").write_text(json.dumps({
        "config_sha256": config["config_sha256"], "mass_GeV": mass, "has_sensitivity": True,
        "u2_min": 1e-8, "u2_max": 1e-4}))
    finalize_and_compact(run_dir, {"name": "central"}, config, [mass])
    assert not vector_dir.exists()
    assert not geometry_dir.exists()
    assert (run_dir / "results" / "mS_0p500.json").exists()
    assert json.loads((run_dir / "complete.json").read_text())["state"] == "compacted"
    assert validate_retained_completion(run_dir, config, [mass]) is not None


def _grid_variation(tmp_path, **extra):
    grid = tmp_path / "central.dat"
    grid.write_text("grid\n")
    return {"name": "central", "axis": "central", "width_scheme": "winkler", "path": grid,
            "sha256": hashlib.sha256(grid.read_bytes()).hexdigest(), **extra}


def _committed_code_hashes(head):
    return {rel: hashlib.sha256(subprocess.check_output(["git", "show", f"{head}:{rel}"], cwd=repo_root())).hexdigest()
            for rel in CODE_INPUTS}


def _embree_available():
    """campaign_config records the pinned Embree backend and refuses to run
    without it (grendel.band.scalar.campaign.ray_backend_provenance)."""
    try:
        import embreex  # noqa: F401
    except ImportError:
        return False
    return True


needs_embree = pytest.mark.skipif(not _embree_available(), reason="the BC4 campaign requires embreex==4.4.0")


def _campaign_code_committed():
    """The retained-config checks hash the campaign code as committed at HEAD,
    so they need a checkout in which every input file is committed."""
    return all(subprocess.run(["git", "cat-file", "-e", f"HEAD:{rel}"], cwd=repo_root(),
                              capture_output=True).returncode == 0 for rel in CODE_INPUTS)


@needs_embree
def test_collector_uses_verified_recorded_producer_config(tmp_path, monkeypatch):
    if not _campaign_code_committed():
        pytest.skip("needs the campaign code committed at HEAD")
    variation = _grid_variation(tmp_path)
    masses = [0.5]
    producer_head = git_head()
    committed = _committed_code_hashes(producer_head)
    monkeypatch.setattr(uncertainty_band, "code_hashes", lambda files: committed)
    recorded = campaign_config(variation, masses, 10, 3, 42)
    run_dir = tmp_path / "runs" / "central"
    run_dir.mkdir(parents=True)
    (run_dir / "run_metadata.json").write_text(json.dumps({**recorded, "producer_git_head_at_start": producer_head}))
    # The working tree may have moved on; the recorded state is checked
    # against the named commit, not against the files now on disk.
    monkeypatch.setattr(uncertainty_band, "code_hashes", lambda files: {"collector.py": "new-state"})
    loaded = load_recorded_config(run_dir, variation, masses, n_pool=10, n_samples=3, seed=42)
    assert loaded["config_sha256"] == recorded["config_sha256"]
    assert loaded["code_sha256"] == recorded["code_sha256"]
    assert loaded["producer_git_head_at_start"] == producer_head


@needs_embree
def test_collector_rejects_tampered_recorded_config(tmp_path):
    if not _campaign_code_committed():
        pytest.skip("needs the campaign code committed at HEAD")
    variation = _grid_variation(tmp_path)
    recorded = campaign_config(variation, [0.5], 10, 3, 42)
    recorded["n_parent_pool"] = 11
    run_dir = tmp_path / "runs" / "central"
    run_dir.mkdir(parents=True)
    (run_dir / "run_metadata.json").write_text(json.dumps({**recorded, "producer_git_head_at_start": git_head()}))
    with pytest.raises(RuntimeError, match="retained config checksum mismatch"):
        load_recorded_config(run_dir, variation, [0.5], n_pool=10, n_samples=3, seed=42)


def _fake_templates(root, seed=7, payload=b"x"):
    root.mkdir(parents=True, exist_ok=True)
    (root / "MANIFEST.json").write_text(json.dumps({
        "decay_model": "exhad:scalar-1809", "exhad_commit": "abc", "exhad_variation": "central",
        "pythia_version": "8.317", "n_templates": 3, "base_seed": seed}))
    (root / "templates_0p500.npz").write_bytes(payload)
    (root / "templates_1p000.npz").write_bytes(payload + b"y")


@needs_embree
def test_template_variation_config_pins_the_bundles(tmp_path):
    templates = tmp_path / "templates"
    _fake_templates(templates)
    analytic = _grid_variation(tmp_path)
    exhad = {**analytic, "name": EXHAD_DECAY_VARIATION, "axis": "decay_model", "width_scheme": EXHAD_DECAY_SCHEME,
             "decay_templates": templates}
    assert config_keys(analytic) == CONFIG_KEYS
    assert config_keys(exhad)[-1] == "decay_templates"
    plain = campaign_config(analytic, [0.5, 1.0], 10, 3, 42)
    templated = campaign_config(exhad, [0.5, 1.0], 10, 3, 42)
    assert "decay_templates" not in plain
    assert templated["decay_templates"]["n_bundles"] == 2
    assert templated["decay_templates"]["base_seed"] == 7
    assert templated["decay_templates"]["decay_model"] == "exhad:scalar-1809"
    assert templated["config_sha256"] != plain["config_sha256"]
    # A changed bundle changes the identity.
    templates_provenance.cache_clear()
    (templates / "templates_1p000.npz").write_bytes(b"different")
    retemplated = campaign_config(exhad, [0.5, 1.0], 10, 3, 42)
    assert retemplated["decay_templates"]["bundle_tree_sha256"] != templated["decay_templates"]["bundle_tree_sha256"]
    assert retemplated["config_sha256"] != templated["config_sha256"]
    templates_provenance.cache_clear()


@needs_embree
def test_template_variation_retained_config_round_trips(tmp_path, monkeypatch):
    if not _campaign_code_committed():
        pytest.skip("needs the campaign code committed at HEAD")
    templates = tmp_path / "templates"
    _fake_templates(templates)
    exhad = _grid_variation(tmp_path, decay_templates=templates)
    exhad.update(name=EXHAD_DECAY_VARIATION, axis="decay_model", width_scheme=EXHAD_DECAY_SCHEME)
    producer_head = git_head()
    monkeypatch.setattr(uncertainty_band, "code_hashes", lambda files: _committed_code_hashes(producer_head))
    templates_provenance.cache_clear()
    recorded = campaign_config(exhad, [0.5], 10, 3, 42)
    run_dir = tmp_path / "runs" / EXHAD_DECAY_VARIATION
    run_dir.mkdir(parents=True)
    (run_dir / "run_metadata.json").write_text(json.dumps({**recorded, "producer_git_head_at_start": producer_head}))
    loaded = load_recorded_config(run_dir, exhad, [0.5], n_pool=10, n_samples=3, seed=42)
    assert loaded["decay_templates"] == recorded["decay_templates"]
    # The same run without its template provenance is incomplete.
    stripped = {k: v for k, v in recorded.items() if k != "decay_templates"}
    (run_dir / "run_metadata.json").write_text(json.dumps({**stripped, "producer_git_head_at_start": producer_head}))
    with pytest.raises(RuntimeError, match="incomplete retained config"):
        load_recorded_config(run_dir, exhad, [0.5], n_pool=10, n_samples=3, seed=42)
    templates_provenance.cache_clear()
