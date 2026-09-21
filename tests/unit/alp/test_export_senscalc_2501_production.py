"""Tests for the staged SensCalc 2501 production export boundary."""

import json
from pathlib import Path

import pytest

from grendel.models.alp.tools import export_senscalc_2501_production as export


def _write_complete_export(path):
    path.mkdir()
    hashes = {}
    for name in export.EXPECTED_OUTPUT_FILES:
        if name == "PRODUCTION_EXPORT_MANIFEST.json":
            continue
        product = path / name
        product.write_text(f"test product: {name}\n")
        hashes[name] = export.sha256(product)
    (path / "PRODUCTION_EXPORT_MANIFEST.json").write_text(json.dumps({
        "senscalc_tag": export.SENSCALC_TAG,
        "senscalc_commit": export.SENSCALC_COMMIT,
        "phenomenology": "arXiv:2501.04525",
        "wolfram_system_id": "test-system",
        "source_sha256": {
            key: expected_hash
            for key, (_, expected_hash) in export.SOURCE_FILES.items()
        },
        "output_sha256": hashes,
    }))


def test_verify_sources_accepts_pinned_checkout():
    checkout = Path("/tmp/SensCalc-review")
    if not checkout.is_dir():
        pytest.skip("pinned SensCalc review checkout is not present")
    verified = export.verify_sources(checkout)
    assert set(verified) == set(export.SOURCE_FILES)


def test_validate_export_accepts_complete_hashed_bundle(tmp_path):
    staged = tmp_path / "staged"
    _write_complete_export(staged)
    export.validate_export(staged)


def test_validate_export_rejects_modified_product(tmp_path):
    staged = tmp_path / "staged"
    _write_complete_export(staged)
    (staged / "fragmentation_probability_coefficients_bnt.csv").write_text(
        "modified\n"
    )
    with pytest.raises(export.InputError, match="hash mismatch"):
        export.validate_export(staged)


def test_wolfram_exporter_maps_association_values():
    source = export.EXPORTER.read_text()
    assert "sourcePaths = Map[" in source
    assert '"source_sha256" -> Map[sha256, sourcePaths]' in source


def test_install_export_replaces_destination_bundle(tmp_path):
    staged = tmp_path / "staged"
    destination = tmp_path / "published"
    _write_complete_export(staged)
    destination.mkdir()
    stale = destination / "fragmentation_probability_coefficients_bnt.csv"
    stale.write_text("stale\n")

    export.install_export(staged, destination)

    output = destination / "fragmentation_probability_coefficients_bnt.csv"
    assert output.read_text().startswith("test product")
    assert not any(staged.iterdir())
