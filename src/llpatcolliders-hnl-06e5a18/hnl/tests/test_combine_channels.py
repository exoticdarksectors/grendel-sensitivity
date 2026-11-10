"""Regression tests for production.combine_channels.

Combination is strict by default: a channel CSV that does not exist (job never
ran) or cannot be parsed fails the combine, while a zero-byte CSV is the
explicit kinematically-closed sentinel and is always accepted. Non-strict mode
is for deliberate partial runs only.
"""

import numpy as np
import pytest

from production import combine_channels


def test_channels_includes_both_tau_variants():
    """tau and induced_tau both in the default CHANNELS list."""
    assert "tau" in combine_channels.CHANNELS
    assert "induced_tau" in combine_channels.CHANNELS
    for ch in ("Bmeson", "Dmeson", "Bc", "Bbaryon", "Kmeson", "WZ"):
        assert ch in combine_channels.CHANNELS


def _write_rows(base, flavor, channel, mass_label, rows):
    ch_dir = base / flavor / channel
    ch_dir.mkdir(parents=True, exist_ok=True)
    np.savetxt(ch_dir / f"mN_{mass_label}.csv", np.asarray(rows),
               delimiter=",", fmt="%.8e")


def test_strict_rejects_missing_channel(tmp_path, monkeypatch):
    """Default (strict): a channel with no CSV at all fails the combine."""
    monkeypatch.setattr(combine_channels, "OUTPUT_BASE", tmp_path)
    _write_rows(tmp_path, "Ue", "Bmeson", "1p000",
                [[1.0, 50.0, 1.0, 2.0, 3.0]])
    with pytest.raises(RuntimeError, match="missing=.*Dmeson"):
        combine_channels.combine_for_point("Ue", 1.0)


def test_strict_rejects_malformed_channel(tmp_path, monkeypatch):
    """Default (strict): a CSV with fewer than 5 columns fails the combine."""
    monkeypatch.setattr(combine_channels, "OUTPUT_BASE", tmp_path)
    _write_rows(tmp_path, "Ue", "Bmeson", "1p000", [[1.0, 50.0, 1.0]])
    with pytest.raises(RuntimeError, match="bad=.*Bmeson"):
        combine_channels.combine_for_point("Ue", 1.0, channels=["Bmeson"])


def test_strict_accepts_zero_byte_sentinels(tmp_path, monkeypatch):
    """Zero-byte CSVs (kinematically-closed channels) pass strict combination."""
    monkeypatch.setattr(combine_channels, "OUTPUT_BASE", tmp_path)
    # Dmeson is closed (D -> tau N closed for Utau on the full grid).
    (tmp_path / "Utau" / "Dmeson").mkdir(parents=True)
    (tmp_path / "Utau" / "Dmeson" / "mN_0p300.csv").write_text("")
    _write_rows(tmp_path, "Utau", "Bmeson", "0p300",
                [[7.0, 80.0, 0.0, 0.0, 70.0]])

    n = combine_channels.combine_for_point(
        "Utau", 0.3, channels=["Dmeson", "Bmeson"])
    assert n == 1


def test_non_strict_skips_missing_channels(tmp_path, monkeypatch):
    """strict=False (partial runs): missing channels contribute zero rows."""
    monkeypatch.setattr(combine_channels, "OUTPUT_BASE", tmp_path)
    _write_rows(tmp_path, "Ue", "Bmeson", "1p000",
                [[1.0, 50.0, 1.0, 2.0, 3.0],
                 [2.0, 60.0, 4.0, 5.0, 6.0]])

    n = combine_channels.combine_for_point("Ue", 1.0, strict=False)
    assert n == 2

    out = np.loadtxt(tmp_path / "Ue" / "combined" / "mN_1p000.csv",
                     delimiter=",")
    if out.ndim == 1:
        out = out.reshape(1, -1)
    assert out.shape == (2, 5)
