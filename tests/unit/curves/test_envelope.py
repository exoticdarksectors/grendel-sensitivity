"""The single-source variation envelope drawn around a BC4/BC10 island."""
import csv
from pathlib import Path

import pytest

from curves.plot.bc_panel import envelope_segments, load_variation_envelope

CONFIG = {"prefix": "u2"}
FIELDS = ["mass_GeV", "has_sensitivity", "envelope_definition", "u2_min_central", "u2_min_envelope_lo",
          "u2_min_envelope_hi", "u2_min_open", "u2_max_central", "u2_max_envelope_lo", "u2_max_envelope_hi",
          "u2_max_open"]


def _write_envelope(path: Path, central_at_first_mass: float = 1.0e-8) -> None:
    rows = [[0.2, True, "single_source_variation_envelope", central_at_first_mass, 0.8e-8, 1.2e-8, False,
             1.0e-4, 0.8e-4, 1.2e-4, False],
            [0.3, False, "single_source_variation_envelope", "", "", "", False, "", "", "", False],
            [0.4, True, "single_source_variation_envelope", 2.0e-8, 1.7e-8, 2.3e-8, False, 2.0e-4, 1.7e-4,
             2.3e-4, False]]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(FIELDS)
        writer.writerows(rows)


def _rewrite(path: Path, edit) -> None:
    rows = list(csv.DictReader(path.open()))
    rows = edit(rows)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def _central_rows():
    return [{"mass_GeV": 0.2, "has_sensitivity": True, "c_min": 1.0e-8, "c_max": 1.0e-4},
            {"mass_GeV": 0.3, "has_sensitivity": False, "c_min": float("nan"), "c_max": float("nan")},
            {"mass_GeV": 0.4, "has_sensitivity": True, "c_min": 2.0e-8, "c_max": 2.0e-4}]


def test_envelope_does_not_bridge_insensitive_rows(tmp_path):
    envelope = tmp_path / "envelope.csv"
    _write_envelope(envelope)
    rows = load_variation_envelope(envelope, CONFIG, _central_rows())
    assert [len(segment) for segment in envelope_segments(rows, "min")] == [1, 1]
    assert [len(segment) for segment in envelope_segments(rows, "max")] == [1, 1]


def test_envelope_rejects_a_different_central_curve(tmp_path):
    envelope = tmp_path / "envelope.csv"
    _write_envelope(envelope, central_at_first_mass=1.1e-8)
    with pytest.raises(ValueError, match="differs from the central curve"):
        load_variation_envelope(envelope, CONFIG, _central_rows())


def test_envelope_accepts_nan_for_an_open_central_edge(tmp_path):
    envelope = tmp_path / "envelope.csv"
    _write_envelope(envelope)

    def open_first_max(rows):
        rows[0].update(u2_max_central="", u2_max_envelope_lo="", u2_max_envelope_hi="", u2_max_open="True")
        return rows

    _rewrite(envelope, open_first_max)
    central = _central_rows()
    central[0]["c_max"] = 1.0e-2
    central[0]["c_max_open"] = True
    loaded = load_variation_envelope(envelope, CONFIG, central)
    assert loaded[0]["max_open"] is True


def test_envelope_rejects_an_open_state_mismatch(tmp_path):
    envelope = tmp_path / "envelope.csv"
    _write_envelope(envelope)

    def open_first_max(rows):
        rows[0]["u2_max_open"] = "True"
        return rows

    _rewrite(envelope, open_first_max)
    with pytest.raises(ValueError, match="open differs from the central"):
        load_variation_envelope(envelope, CONFIG, _central_rows())


def test_sparse_envelope_breaks_at_an_intervening_gap(tmp_path):
    envelope = tmp_path / "envelope.csv"
    _write_envelope(envelope)
    _rewrite(envelope, lambda rows: [rows[0], rows[2]])
    loaded = load_variation_envelope(envelope, CONFIG, _central_rows())
    assert [len(segment) for segment in envelope_segments(loaded, "min")] == [1, 1]
    assert [len(segment) for segment in envelope_segments(loaded, "max")] == [1, 1]


@pytest.mark.parametrize("diagnostic", ("variation_open", "variation_missing", "numerically_unresolved"))
def test_envelope_splits_at_unsupported_campaign_diagnostics(diagnostic):
    def row(mass):
        return {"mass_GeV": mass, "has_sensitivity": True, "min_open": False, "min_variation_open": False,
                "min_variation_missing": False, "min_numerically_unresolved": False, "min_lo": 0.8e-8,
                "min_hi": 1.2e-8}

    rows = [row(0.2), row(0.3), row(0.4)]
    rows[1][f"min_{diagnostic}"] = True
    assert [len(segment) for segment in envelope_segments(rows, "min")] == [1, 1]
