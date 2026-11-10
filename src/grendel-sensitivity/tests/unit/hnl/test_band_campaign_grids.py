"""The HNL FONLL campaign verifies its grids and requires the complete campaign."""
import json

import pytest

from grendel.band.fonll_variations import grid_path
from grendel.band.hnl.campaign import discover_variations
from grendel.io.atomic import sha256_file

SCALES = ["scale_muR2_muF2", "scale_muR0p5_muF0p5", "scale_muR2_muF1", "scale_muR1_muF2",
          "scale_muR0p5_muF1", "scale_muR1_muF0p5"]
MASSES = {"bottom": ["mass_4p5", "mass_5"], "charm": ["mass_1p3", "mass_1p7"]}


def _campaign(grid_dir, drop=()):
    grid_dir.mkdir()
    grids = []
    for quark in ("bottom", "charm"):
        tags = (["central"] + SCALES + [f"pdf_{i:04d}" for i in range(1, 101)] + MASSES[quark])
        for tag in tags:
            if (quark, tag) in drop:
                continue
            path = grid_path(grid_dir, tag, quark)
            path.write_text(f"{quark} {tag}\n")
            kind = tag.split("_")[0] if tag != "central" else "central"
            grids.append({"quark": quark, "variation_tag": tag, "variation_kind": kind,
                          "sha256": sha256_file(path)})
    (grid_dir / "variation_manifest.json").write_text(json.dumps({"grids": grids}))
    return grid_dir


def test_complete_campaign_is_discovered_and_selected(tmp_path):
    grid_dir = _campaign(tmp_path / "grids")
    variations = discover_variations(grid_dir)
    axes = [v["axis"] for v in variations]
    assert {axis: axes.count(axis) for axis in set(axes)} == {"central": 1, "scale": 6, "pdf": 100, "mass": 4}
    assert [v["name"] for v in discover_variations(grid_dir, ["central", "mb_up"])] == ["central", "mb_up"]


def test_modified_grid_fails_its_checksum(tmp_path):
    grid_dir = _campaign(tmp_path / "grids")
    grid_path(grid_dir, "pdf_0042", "charm").write_text("edited\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        discover_variations(grid_dir, ["central"])


def test_missing_grid_fails_loudly(tmp_path):
    grid_dir = _campaign(tmp_path / "grids", drop={("bottom", "mass_5")})
    with pytest.raises(FileNotFoundError, match="mass_5"):
        discover_variations(grid_dir)
    grid_dir = _campaign(tmp_path / "grids2")
    grid_path(grid_dir, "scale_muR2_muF1", "bottom").unlink()
    with pytest.raises(FileNotFoundError):
        discover_variations(grid_dir)


def test_partial_replica_set_is_refused(tmp_path):
    grid_dir = _campaign(tmp_path / "grids", drop={("bottom", "pdf_0100"), ("charm", "pdf_0100")})
    with pytest.raises(ValueError, match="incomplete FONLL campaign"):
        discover_variations(grid_dir)
