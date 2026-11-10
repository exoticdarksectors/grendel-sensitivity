"""FONLL build and campaign plumbing without FONLL: the Linux make overrides,
the reuse of raw quark grids, the serial charm calibration and the patch pins
recorded in the manifest."""
import importlib
import platform

import pytest

from grendel.production.fonll_grids import generate as gen, install

# the macOS build the grids were produced with: its compiler line, no trapfpe.o
LINUX_OVERRIDES = ("F77=gfortran -Wall -fno-automatic -fallow-argument-mismatch", "SYSOBJ=")


@pytest.mark.parametrize("system, expected", [("Linux", True), ("Darwin", False)])
def test_make_matches_the_macos_build_on_linux_only(monkeypatch, system, expected):
    monkeypatch.setattr(platform, "system", lambda: system)
    try:
        make = importlib.reload(install).MAKE
        assert [override in make for override in LINUX_OVERRIDES] == [expected, expected]
    finally:
        monkeypatch.undo()
        importlib.reload(install)


def _generate(tmp_path, monkeypatch, raw_rows):
    """generate() with --reuse-existing-grids over a raw grid of ``raw_rows`` rows."""
    ws = gen.Workspace(fonll=tmp_path / "fonll", out=tmp_path / "out", run=tmp_path / "run", logs=tmp_path / "logs")
    grid_file = ws.run / "nlo_bottom_central" / "nl_b_.out"
    grid_file.parent.mkdir(parents=True)
    grid_file.write_text("header\n" + "row\n" * raw_rows)
    built = []

    def build(ws_, run_dir, grid, *args, **kwargs):
        built.append(grid.exists())          # an incomplete grid is removed first (fonllgridlha appends)
        grid.write_text("header\n" + "row\n" * (len(gen.Y_VALUES) * gen.FONLL_GRID_NPT))

    monkeypatch.setattr(gen, "build_quark_grid", build)
    monkeypatch.setattr(gen, "run_fragmentation", lambda ws_, pdf, quark, run_dir, grid, mode, param, points, label:
                        ({gen.point_key(pt, y): 1.0 for pt, y in points}, run_dir / "frag.dat"))
    summary = gen.generate(ws, "nlo", "bottom", gen.central_variation("bottom"), reuse_existing_grids=True,
                           grid_workers=1)
    return built, summary


def test_complete_raw_grid_is_reused(tmp_path, monkeypatch):
    built, summary = _generate(tmp_path, monkeypatch, len(gen.Y_VALUES) * gen.FONLL_GRID_NPT)
    assert built == [] and summary["raw_grid_reused"]


def test_raw_grid_cut_on_a_rapidity_row_is_rebuilt(tmp_path, monkeypatch):
    built, summary = _generate(tmp_path, monkeypatch, (len(gen.Y_VALUES) - 1) * gen.FONLL_GRID_NPT)
    assert built == [False] and not summary["raw_grid_reused"]


def test_calibrate_only_fits_the_charm_weight_and_exits(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(gen, "load_or_fit_charm_feeddown_weight", lambda ws: calls.append(ws.out))
    monkeypatch.setattr(gen, "generate", lambda *args, **kwargs: pytest.fail("no grid is generated"))
    assert gen.main(["--calibrate-only", "--out-dir", str(tmp_path / "out")]) == 0
    assert calls == [tmp_path / "out"]


def test_manifest_records_the_pinned_patch_hashes():
    for rel, pin in install.PINS["files"].items():
        assert f"{rel}: patched sha256 {pin['sha256_patched']}" in gen.PATCH_NOTES
