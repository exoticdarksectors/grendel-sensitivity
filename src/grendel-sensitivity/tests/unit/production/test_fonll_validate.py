"""The FONLL closure checks without FONLL: the gate compares every column of
the regenerated central, including its pT = 0 cells, and the public bottom
point is fragmented at (5 GeV, 0), which is not a node of the meson grid."""
import pytest

from grendel.production.fonll_grids import generate as gen, validate

REFERENCE = validate.REFERENCE_DIR / gen.grid_filename("nlo", "central", "bottom")


def _workspace(tmp_path):
    return gen.Workspace(fonll=tmp_path / "fonll", out=tmp_path / "out", run=tmp_path / "run", logs=tmp_path / "logs")


def _gate(tmp_path, monkeypatch, row=None, edit=None):
    """gate_quark on a copy of the tracked central with data row ``row`` edited."""
    lines = REFERENCE.read_text().splitlines()
    if row is not None:
        data = [i for i, line in enumerate(lines) if not line.startswith("#")]
        lines[data[row]] = " ".join(edit(lines[data[row]].split()))
    fresh = tmp_path / "fresh.dat"
    fresh.write_text("\n".join(lines) + "\n")
    monkeypatch.setattr(gen, "generate", lambda *args, **kwargs: {"path": str(fresh)})
    return validate.gate_quark(_workspace(tmp_path), "bottom", validate.REFERENCE_DIR, 1e-3, 1)


def test_gate_passes_the_tracked_central(tmp_path, monkeypatch):
    assert _gate(tmp_path, monkeypatch)


def test_gate_fails_garbage_where_the_reference_is_zero(tmp_path, monkeypatch):
    assert not _gate(tmp_path, monkeypatch, row=0, edit=lambda c: [c[0], c[1], "1.0e+05"])   # pT = 0, y = -3


def test_gate_fails_moved_nodes(tmp_path, monkeypatch):
    assert not _gate(tmp_path, monkeypatch, row=150, edit=lambda c: [c[0], "0.01", c[2]])


def test_public_bottom_point_is_fragmented_at_the_requested_point(tmp_path, monkeypatch):
    ws = _workspace(tmp_path)
    with pytest.raises(FileNotFoundError, match="generate --pdf cteq66 --quark bottom"):
        validate.public_points(ws, 0.05)
    grid = ws.run / "cteq66_bottom_central" / "ct_b_.out"
    grid.parent.mkdir(parents=True)
    grid.write_text("raw grid\n")
    calls = []

    def fragment(ws_, pdf_key, quark, run_dir, grid_file, frag_mode, frag_param, points, label):
        calls.append((pdf_key, quark, grid_file, frag_mode, frag_param, points))
        return {gen.point_key(5.0, 0.0): validate.PUBLIC_BOTTOM_PT5_Y0 * 1.001}, run_dir / "frag.dat"

    monkeypatch.setattr(gen, "run_fragmentation", fragment)
    monkeypatch.setattr(gen, "fit_charm_d0_feeddown_weight", lambda ws_: {
        "max_abs_relative_difference": 0.002, "vector_to_direct_weight": 1.2, "points": []})
    assert validate.public_points(ws, 0.05)
    assert calls == [("cteq66", "bottom", grid, 2, 24.2, [(5.0, 0.0)])]
