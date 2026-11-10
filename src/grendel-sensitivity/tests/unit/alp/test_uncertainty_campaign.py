"""The BC10 variation campaign: the variation list, seeds, the single-source
envelope, topology bookkeeping and the compaction of completed runs."""
import json

import numpy as np
import pandas as pd
import pytest

from grendel.band.alp import variations as campaign_definitions
from grendel.band.alp.campaign import _compact_completed_run, _scratch_relative, headline, load_mass_grid
from grendel.band.alp.combine import combine_band
from grendel.band.alp.variations import all_variations, numerical_control_variations, stable_seed
from grendel.io.atomic import sha256_file

# Every axis of the campaign, one variation each unless noted.
DEFINITIONS = ([("central", "central")] + [(f"scale_{i}", "scale") for i in range(6)]
               + [(f"pdf_{i}", "pdf") for i in range(100)] + [("mb_dn", "mb"), ("mb_up", "mb")]
               + [(f"gg_{q}", "decay_gg") for q in "uds"] + [("cbs_down", "cbs"), ("cbs_up", "cbs")]
               + [("central_repeat_1", "numerical_control"), ("central_repeat_2", "numerical_control")])


def test_stable_seed_is_reproducible_and_member_specific():
    assert stable_seed("pdf_0001") == stable_seed("pdf_0001")
    assert stable_seed("pdf_0001") != stable_seed("pdf_0002")
    assert 0 < stable_seed("scale_muR2_muF2") < 2**32


def test_numerical_controls_have_independent_production_and_reco_seeds(tmp_path):
    grid = tmp_path / "central.dat"
    grid.write_text("grid")
    controls = numerical_control_variations(grid)
    assert [item["name"] for item in controls] == ["central_repeat_1", "central_repeat_2"]
    assert len({item["production_seed"] for item in controls}) == 2
    assert len({item["reco_seed_offset"] for item in controls}) == 2
    assert all(item["production_seed"] != item["reco_seed_offset"] for item in controls)


def test_campaign_enumerates_every_axis(tmp_path, monkeypatch):
    grid = tmp_path / "central.dat"
    grid.write_text("grid")
    fonll = [{"name": "central", "axis": "central", "campaign_axis": "fonll", "grid_path": str(grid)}]
    for axis, count in (("scale", 6), ("pdf", 100), ("mb", 2)):
        fonll.extend({"name": f"{axis}_{index}", "axis": axis, "campaign_axis": "fonll", "grid_path": str(grid)}
                     for index in range(count))
    monkeypatch.setattr(campaign_definitions, "discover_fonll_variations", lambda _: fonll)
    variations = all_variations(tmp_path)
    counts = {}
    for item in variations:
        counts[item["axis"]] = counts.get(item["axis"], 0) + 1
    assert counts == {"central": 1, "scale": 6, "pdf": 100, "mb": 2, "decay_gg": 3, "cbs": 2, "numerical_control": 2}
    assert len({item["name"] for item in variations}) == len(variations)


def test_mass_grid_file_is_hashed_and_requires_strict_order(tmp_path):
    path = tmp_path / "dense_grid.csv"
    pd.DataFrame({"mass_GeV": [1.18, 1.19, 1.20]}).to_csv(path, index=False)
    result = load_mass_grid(path)
    assert result["n_masses"] == 3
    assert result["masses_GeV"] == [1.18, 1.19, 1.20]
    assert result["source_path"] == str(path.resolve())
    assert result["source_sha256"] == sha256_file(path)
    assert len(result["canonical_sha256"]) == 64
    pd.DataFrame({"mass_GeV": [1.18, 1.20, 1.19]}).to_csv(path, index=False)
    with pytest.raises(ValueError, match="strictly increasing"):
        load_mass_grid(path)


def test_registry_paths_are_scratch_relative(tmp_path):
    scratch = tmp_path / "scratch"
    artifact = scratch / "runs" / "central" / "variation.complete.json"
    artifact.parent.mkdir(parents=True)
    artifact.write_text("{}")
    assert _scratch_relative(artifact, scratch) == "runs/central/variation.complete.json"
    with pytest.raises(ValueError, match="outside the campaign directory"):
        _scratch_relative(tmp_path / "elsewhere", scratch)


def _rows(mass, sensitive_names=None, shifts=None, xc=-8.0):
    rows = []
    for name, axis in DEFINITIONS:
        sensitive = True if sensitive_names is None else name in sensitive_names
        shift = (shifts or {}).get(name, 0.0)
        rows.append({"mass_GeV": mass, "variation": name, "axis": axis, "has_sensitivity": sensitive,
                     "invf_min": 10.0 ** (xc + shift) if sensitive else np.nan,
                     "invf_max": 10.0 ** (xc + 1.0 + shift) if sensitive else np.nan,
                     "invf_min_open": False, "invf_max_open": False})
    return pd.DataFrame(rows)


def test_exact_campaign_builds_single_source_variation_envelope():
    xc = -8.0
    pdf_shifts = np.r_[np.linspace(-0.02, 0.02, 99), 1.0]
    shifts = {f"scale_{i}": s for i, s in enumerate((-0.10, -0.05, -0.01, 0.01, 0.05, 0.10))}
    shifts.update({f"pdf_{i}": s for i, s in enumerate(pdf_shifts)})
    shifts.update({"mb_dn": -0.03, "mb_up": 0.03, "gg_u": -0.04, "gg_d": 0.0, "gg_s": 0.04,
                   "cbs_down": -0.05, "cbs_up": 0.05, "central_repeat_1": 0.12, "central_repeat_2": -0.007})
    rows = _rows(1.0, shifts=shifts, xc=xc)
    band = combine_band(rows).iloc[0]
    p16, p84 = np.quantile(pdf_shifts, [0.16, 0.84])
    assert band["envelope_definition"] == "single_source_variation_envelope"
    assert np.isclose(band["invf_min_pdf_p16"], 10.0 ** (xc + p16))
    assert np.isclose(band["invf_min_pdf_p84"], 10.0 ** (xc + p84))
    assert np.isclose(band["invf_min_pdf_log10_std"], np.std(pdf_shifts, ddof=1))
    assert np.isclose(band["invf_min_envelope_lo"], 10.0 ** (xc - 0.10))
    assert np.isclose(band["invf_min_envelope_hi"], 10.0 ** (xc + 0.10))
    assert band["invf_min_envelope_lo_source"] == "scale"
    assert band["invf_min_envelope_hi_source"] == "scale"
    assert band["invf_min_pdf_n_finite"] == 100
    assert np.isclose(band["invf_min_repeat_max_abs_dex"], 0.12)
    assert bool(band["invf_min_repeat_not_subdominant"])
    assert not bool(band["invf_min_variation_missing"])

    reference = pd.DataFrame([{"mass_GeV": 1.0, "has_sensitivity": True, "invf_min": 2.0 * 10.0**xc,
                               "invf_max": 2.0 * 10.0**(xc + 1.0), "invf_min_open": False, "invf_max_open": False}])
    rebased = combine_band(rows, reference).iloc[0]
    assert rebased["invf_min_campaign_central"] == pytest.approx(10.0**xc)
    assert rebased["invf_min_central"] == pytest.approx(2.0 * 10.0**xc)
    assert rebased["invf_min_envelope_lo"] == pytest.approx(2.0 * 10.0**(xc - 0.10))
    assert rebased["invf_min_envelope_hi"] == pytest.approx(2.0 * 10.0**(xc + 0.10))
    assert rebased["envelope_reference"] == "canonical_high_statistics_central"
    assert bool(rebased["canonical_rebase_topology_compatible"])


def test_missing_variation_boundary_is_explicit():
    names = {name for name, _ in DEFINITIONS} - {"scale_0"}
    band = combine_band(_rows(1.0, sensitive_names=names)).iloc[0]
    assert bool(band["invf_min_variation_missing"])
    assert bool(band["invf_max_variation_missing"])
    assert np.isfinite(band["invf_min_envelope_lo"])


def test_central_insensitive_row_remains_an_explicit_gap():
    band = combine_band(_rows(0.54, sensitive_names=set())).iloc[0]
    assert not bool(band["has_sensitivity"])
    assert not bool(band["any_variation_sensitive"])
    assert not bool(band["any_halo_variation_sensitive"])
    assert np.isnan(band["invf_min_envelope_lo"])
    assert np.isnan(band["invf_max_envelope_hi"])


def test_numerical_repeat_topology_difference_is_preserved_at_central_gap():
    band_frame = combine_band(_rows(3.30, sensitive_names={"central_repeat_2"}))
    band = band_frame.iloc[0]
    assert not bool(band["has_sensitivity"])
    assert not bool(band["any_halo_variation_sensitive"])
    assert bool(band["numerical_control_any_sensitive"])
    assert not bool(band["numerical_control_all_sensitive"])
    assert band["numerical_control_n_sensitive"] == 1
    assert band["numerical_control_sensitive_variations"] == "central_repeat_2"
    assert bool(band["numerical_control_topology_differs"])
    assert band["numerical_control_topology_difference_variations"] == "central_repeat_2"
    assert np.isnan(band["invf_min_central_repeat_1"])
    assert band["invf_min_central_repeat_2"] == pytest.approx(1e-8)
    assert not bool(band["numerical_control_included_in_halo"])
    assert headline(band_frame)["numerical_control_topology"] == {
        "n_differences": 1, "difference_masses_GeV": [3.3], "difference_variations_by_mass": {"3.3": "central_repeat_2"}}


def test_physical_variation_topology_difference_is_named_at_central_gap():
    band_frame = combine_band(_rows(3.30, sensitive_names={"gg_u"}))
    band = band_frame.iloc[0]
    assert not bool(band["has_sensitivity"])
    assert bool(band["any_halo_variation_sensitive"])
    assert bool(band["halo_restores_sensitivity"])
    assert band["halo_restores_sensitivity_variations"] == "gg_u"
    assert not bool(band["halo_removes_sensitivity"])
    assert bool(band["halo_topology_differs"])
    assert band["halo_topology_difference_variations"] == "gg_u"
    assert np.isnan(band["invf_min_envelope_lo"])
    assert headline(band_frame)["physical_variation_topology"] == {
        "n_differences": 1, "difference_masses_GeV": [3.3], "difference_variations_by_mass": {"3.3": "gg_u"},
        "restored_masses_GeV": [3.3], "removed_masses_GeV": []}


def test_post_validation_compaction_is_atomic_and_idempotent(tmp_path):
    run_dir = tmp_path / "runs" / "pdf_0001"
    vectors = run_dir / "llp_4vectors"
    geometry = run_dir / "analysis" / "geometry_cache"
    vectors.mkdir(parents=True)
    geometry.mkdir(parents=True)
    (vectors / "mA_1.csv").write_bytes(b"vectors")
    (geometry / "geom_1.npz").write_bytes(b"geometry")
    marker = {"variation": {"name": "pdf_0001", "production_mode": "fresh_600k"},
              "sensitivity_csv": str(run_dir / "analysis" / "sensitivity.csv"),
              "production_marker": str(run_dir / "production.complete.json"), "storage_state": "full"}
    completion = run_dir / "variation.complete.json"
    completion.write_text(json.dumps(marker))
    result = _compact_completed_run(completion, marker, run_dir, False)
    assert result["storage_state"] == "compacted"
    assert result["compaction"]["reclaimed_bytes"] == len(b"vectorsgeometry")
    assert not vectors.exists()
    assert not geometry.exists()
    assert json.loads(completion.read_text())["storage_state"] == "compacted"
    assert _compact_completed_run(completion, result, run_dir, False) == result
