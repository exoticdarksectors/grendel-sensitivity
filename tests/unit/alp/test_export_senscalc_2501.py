"""Tests for the staged SensCalc 2501 export boundary."""

import csv
import json

import pytest

from grendel.models.alp.tools import export_senscalc_2501 as export


def _write_complete_export(path, decay_model=export.DEFAULT_DECAY_MODEL):
    path.mkdir()
    width_ids = [
        "mass_GeV",
        "width_coefficient_001",
        "width_coefficient_002",
        "width_coefficient_003",
        "width_coefficient_004",
        "width_coefficient_005",
        "width_coefficient_006",
        "width_coefficient_007",
    ]
    canonical_widths = [
        None,
        "ee",
        "mumu",
        "tautau",
        "gammagamma",
        "nonhadronic_total",
        "hadronic_total",
        "total",
    ]
    raw_rows = [
        [0.5, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0],
        [1.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0],
    ]

    def write_csv(name, header, rows):
        with (path / name).open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(rows)

    write_csv("widths_raw_senscalc.csv", width_ids, raw_rows)
    write_csv(
        "widths_bnt.csv",
        width_ids,
        [[row[0], *(value / 4.0 for value in row[1:])] for row in raw_rows],
    )
    (path / "widths_metadata.json").write_text(json.dumps({
        "columns": [
            {"id": column_id, "canonical_name": canonical_name}
            for column_id, canonical_name in zip(width_ids, canonical_widths)
        ],
        "conversion_raw_to_bnt": 0.25,
    }))

    channel_ids = [f"channel_{index:03d}" for index in range(1, 5)]
    canonical_channels = ["ee", "mumu", "tautau", "gammagamma"]
    write_csv(
        "branching_ratios.csv",
        ["mass_GeV", *channel_ids],
        [[0.5, 0.25, 0.25, 0.25, 0.25], [1.0, 0.2, 0.3, 0.1, 0.4]],
    )
    (path / "decay_channels.json").write_text(json.dumps([
        {"id": channel_id, "canonical_name": canonical_name}
        for channel_id, canonical_name in zip(
            channel_ids, canonical_channels
        )
    ]))
    (path / "matrix_elements.json").write_text(json.dumps([
        {"id": "matrix_element_001", "canonical_name": None}
    ]))

    hashes = {}
    for name in export.EXPECTED_OUTPUT_FILES:
        if name == "EXPORT_MANIFEST.json":
            continue
        product = path / name
        assert product.is_file()
        hashes[name] = export.sha256(product)
    (path / "EXPORT_MANIFEST.json").write_text(json.dumps({
        "senscalc_tag": export.SENSCALC_TAG,
        "senscalc_commit": export.SENSCALC_COMMIT,
        "phenomenology": export.DECAY_MODEL_SPECS[decay_model]["phenomenology"],
        "decay_model": decay_model,
        "wolfram_system_id": "test-system",
        "source_sha256": {
            key: expected_hash
            for key, (_, expected_hash) in export.source_files(decay_model).items()
        },
        "output_sha256": hashes,
    }))


def _refresh_manifest_hash(path, name):
    manifest_path = path / "EXPORT_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["output_sha256"][name] = export.sha256(path / name)
    manifest_path.write_text(json.dumps(manifest))


def test_validate_export_accepts_complete_hashed_bundle(tmp_path):
    staged = tmp_path / "staged"
    _write_complete_export(staged)
    export.validate_export(staged)


def test_validate_export_accepts_exact_2310_structural_schema(tmp_path):
    staged = tmp_path / "staged"
    _write_complete_export(staged, "2310_structural")
    export.validate_export(staged, "2310_structural")


def test_2310_source_files_and_hashes_are_pinned():
    sources = export.source_files("2310_structural")
    assert sources["widths"][0].name.endswith("2310.03524.m")
    assert sources["branching_ratios"][1] == (
        "130aaa7e95c50a3fea8f43a537fcddfe2bf2d4f676a6110bd3f88097da2826cd"
    )
    assert sources["matrix_elements"][1] == (
        "89ef54601fc995418e9967d8269862e27cb86852a847bff1669e8948aa12dfca"
    )


def test_validate_export_rejects_modified_product(tmp_path):
    staged = tmp_path / "staged"
    _write_complete_export(staged)
    (staged / "widths_bnt.csv").write_text("modified\n")
    with pytest.raises(export.InputError, match="hash mismatch"):
        export.validate_export(staged)


def test_validate_export_rejects_wrong_bnt_conversion(tmp_path):
    staged = tmp_path / "staged"
    _write_complete_export(staged)
    bnt_path = staged / "widths_bnt.csv"
    text = bnt_path.read_text().replace("1.0,2.0", "1.0,2.5", 1)
    bnt_path.write_text(text)
    _refresh_manifest_hash(staged, "widths_bnt.csv")
    with pytest.raises(export.InputError, match="not raw/4"):
        export.validate_export(staged)


def test_validate_export_rejects_missing_canonical_width(tmp_path):
    staged = tmp_path / "staged"
    _write_complete_export(staged)
    metadata_path = staged / "widths_metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["columns"][-1]["canonical_name"] = None
    metadata_path.write_text(json.dumps(metadata))
    _refresh_manifest_hash(staged, "widths_metadata.json")
    with pytest.raises(export.InputError, match="canonical width"):
        export.validate_export(staged)


def test_mx_system_id_reads_binary_header(tmp_path):
    source = tmp_path / "input.mx"
    source.write_bytes(
        b"(*This is a Wolfram Language binary dump file.*)\x00"
        b"\x00\x01Windows-x86-64\x00payload"
    )
    assert export.mx_system_id(source) == "Windows-x86-64"


def test_mx_system_id_ignores_plain_text(tmp_path):
    source = tmp_path / "input.txt"
    source.write_text("mass coefficient\n")
    assert export.mx_system_id(source) is None


def test_wolfram_exporter_maps_association_values():
    source = export.EXPORTER.read_text()
    assert "sourcePaths = Map[" in source
    assert '"source_sha256" -> Map[sha256, sourcePaths]' in source


def test_install_export_replaces_destination_bundle(tmp_path):
    staged = tmp_path / "staged"
    destination = tmp_path / "published"
    _write_complete_export(staged)
    destination.mkdir()
    (destination / "widths_bnt.csv").write_text("stale\n")

    export.install_export(staged, destination)

    assert (destination / "widths_bnt.csv").read_text().startswith("mass_GeV")
    assert not any(staged.iterdir())
