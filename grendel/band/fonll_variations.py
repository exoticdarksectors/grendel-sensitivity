"""The coherent FONLL grid campaign behind every production band.

The heavy-flavour grids are produced once, by the FONLL generator in
``grendel.production.fonll_grids``, as a campaign: one central grid, the
six non-central members of the seven-point scale set, the 100 NNPDF4.0
Monte-Carlo replicas, and two heavy-quark-mass grids per quark, all listed
with checksums in ``variation_manifest.json``. Every model's campaign reads
that manifest, verifies the grids it uses, and names its variations from
the entries; the expected counts are asserted so a partial campaign cannot
silently produce a narrower band.
"""
from __future__ import annotations

import json
from pathlib import Path

from ..io.atomic import sha256_file

GRID_STEM = "fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3"
KINDS = ("central", "scale", "pdf", "mass")


def load_manifest(grid_dir) -> dict:
    return json.loads((Path(grid_dir) / "variation_manifest.json").read_text())


def manifest_index(manifest: dict) -> dict[tuple[str, str], dict]:
    """Entries by ``(quark, variation_tag)``."""
    return {(e["quark"], e["variation_tag"]): e for e in manifest["grids"]}


def grid_path(grid_dir, tag: str, quark: str) -> Path:
    return Path(grid_dir) / f"{GRID_STEM}_{tag}_{quark}.dat"


def verified_sha(path: Path, entry: dict, validate_hashes: bool = True) -> str:
    if not path.exists():
        raise FileNotFoundError(path)
    actual = sha256_file(path)
    if validate_hashes and actual != entry["sha256"]:
        raise ValueError(f"checksum mismatch for {path.name}: {actual} != {entry['sha256']}")
    return actual


def bottom_records(grid_dir, *, validate_hashes: bool = True,
                   mass_naming: str = "direction") -> list[dict]:
    """The bottom-grid ensemble as records ``{name, axis, tag, kind, path,
    sha256, muR, muF, lhapdf_member, heavy_quark_mass_GeV,
    trapezoid_integral_pb}``.

    PDF member 0 (the replica mean) is skipped. Mass grids are named by
    direction relative to the central mass (``mb_dn``/``mb_up``) or by their
    manifest tag (``mass_naming="tag"``).
    """
    grid_dir = Path(grid_dir)
    manifest = load_manifest(grid_dir)
    bottom = [g for g in manifest["grids"] if g["quark"] == "bottom"]
    central_entry = next((e for e in bottom if e["variation_tag"] == "central"), None)
    if central_entry is None:
        raise ValueError("FONLL variation manifest has no central bottom grid")
    central_mass = float(central_entry["heavy_quark_mass_GeV"])
    records = []
    for entry in bottom:
        tag, kind = entry["variation_tag"], entry["variation_kind"]
        if tag == "central":
            name, axis = "central", "central"
        elif kind == "scale":
            name, axis = tag, "scale"
        elif kind == "pdf":
            if int(entry.get("lhapdf_member", 0)) == 0:
                continue
            name, axis = tag, "pdf"
        elif kind == "mass":
            if mass_naming == "tag":
                name = tag
            else:
                name = "mb_dn" if float(entry["heavy_quark_mass_GeV"]) < central_mass else "mb_up"
            axis = "mass"
        else:
            continue
        path = grid_path(grid_dir, tag, "bottom")
        records.append({
            "name": name, "axis": axis, "tag": tag, "kind": kind, "path": path,
            "sha256": verified_sha(path, entry, validate_hashes),
            "trapezoid_integral_pb": float(entry["trapezoid_integral_pb"]),
            "muR": entry.get("muR"), "muF": entry.get("muF"),
            "lhapdf_member": entry.get("lhapdf_member"),
            "heavy_quark_mass_GeV": entry.get("heavy_quark_mass_GeV"),
        })
    order = {"central": 0, "scale": 1, "pdf": 2, "mass": 3}
    records.sort(key=lambda r: (order[r["axis"]], r["name"]))
    return records


def assert_counts(records, expected: dict) -> None:
    counts = {axis: sum(r["axis"] == axis for r in records) for axis in expected}
    if counts != expected:
        raise ValueError(f"incomplete FONLL campaign: found {counts}, expected {expected}")


def select_variations(variations, selected):
    """Filter records by name or axis; None keeps everything."""
    if not selected:
        return variations
    wanted = set(selected)
    picked = [r for r in variations if r["name"] in wanted or r["axis"] in wanted]
    if not picked:
        raise ValueError(f"variation selection matched nothing: {sorted(wanted)}")
    return picked
