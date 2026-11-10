"""The FONLL grid tooling without FONLL: enumeration, input cards, the
charm feeddown combination, the manifest partition and the envelope
combiner on synthetic grids, and the patch installer's hash checks."""
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
import pytest

from grendel.production.fonll_grids import combine, generate as gen, install, snapshot

PATCHES = Path(install.__file__).parent / "patches"


# --------------------------------------------------------------- generate --

def test_full_campaign_enumerates_218_grids():
    members = list(range(0, 101))
    per_quark = {q: gen.enumerate_variations(q, members, True, True) for q in gen.QUARKS}
    assert {q: len(v) for q, v in per_quark.items()} == {"bottom": 109, "charm": 109}
    for variations in per_quark.values():
        assert variations[0].kind == "central" and variations[0].tag == "central"
        assert len({v.tag for v in variations}) == len(variations)
        assert sum(v.kind == "scale" for v in variations) == 6      # (1,1) is the central
        assert sum(v.kind == "pdf" for v in variations) == 100      # member 0 is the central
        assert sum(v.kind == "mass" for v in variations) == 2


def test_variation_tags_are_filesystem_safe():
    tags = [v.tag for v in gen.enumerate_variations("bottom", [0, 7], True, True)]
    assert "scale_muR2_muF1" in tags and "scale_muR0p5_muF1" in tags and "pdf_0007" in tags and "mass_4p5" in tags
    assert all(set(t) <= set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_") for t in tags)


def test_parse_pdf_members_always_includes_the_central():
    assert gen.parse_pdf_members("1-3") == [0, 1, 2, 3]
    assert gen.parse_pdf_members("5,1,5") == [0, 1, 5]
    assert gen.parse_pdf_members("") == [0]


def test_grid_input_card_has_the_fonll_layout():
    card = gen.grid_input("nl_b_", 331700, 4.75, ffact=2.0, fren=0.5, y_values=[-1.0, 0.0, 1.0]).splitlines()
    assert card[0] == "nl_b_"
    assert card[1] == card[2] == " 1 7000. 0 0 331700"
    assert card[3] == " 4.75" and card[5] == " 2 0.5"          # (ffact = muF, fren = muR)
    assert card[6:9] == [" -1", " 0", " 1"]
    assert card[-2:] == [" 200 80", " 1"]


def test_frag_input_card_lists_points_and_terminates():
    card = gen.frag_input_points("nl_b_.out", 2, 24.2, [(5.0, 0.0), (10.0, 1.5)]).splitlines()
    assert card[3] == "nl_b_.out" and card[6] == "2" and card[9] == "24.2"
    assert card[-3:] == ["5 0", "10 1.5", "-1 0"]


def test_charm_feeddown_reduces_to_the_direct_stream_at_zero_weight():
    points = [(5.0, 0.0), (20.0, 1.0)]
    direct = {gen.point_key(pt, y): 1.0 + pt for pt, y in points}
    vector = {gen.point_key(pt / float(c["z_collinear"]), y): 100.0 for pt, y in points
              for c in gen.charm_d0_feeddown_channels()}
    combined = gen.combine_charm_d0_feeddown(direct, vector, 0.0, points)
    assert combined == {gen.point_key(pt, y): 1.0 + pt for pt, y in points}
    # With a weight the result is the BR-weighted average of the two streams.
    channels = gen.charm_d0_feeddown_channels()
    br_sum = sum(float(c["branching_fraction"]) for c in channels)
    feed = sum(float(c["branching_fraction"]) * 100.0 / float(c["z_collinear"]) for c in channels)
    w = 1.2
    combined_w = gen.combine_charm_d0_feeddown(direct, vector, w, points)
    assert combined_w[gen.point_key(5.0, 0.0)] == pytest.approx((6.0 + w * feed) / (1 + w * br_sum))


def test_collinear_energy_fractions_are_below_one():
    for c in gen.charm_d0_feeddown_channels():
        assert 0.9 < float(c["z_collinear"]) < 1.0


# ---------------------------------------------------------------- combine --

def _write_grid(path: Path, quark: str, tag: str, member: int, muR: float, muF: float, mass: float,
                scale: float) -> dict:
    ref_pt, ref_y = combine.reference_grid_coords()
    values = scale * (1.0 + ref_pt / 50.0) * np.exp(-ref_y ** 2)
    with path.open("w") as fh:
        fh.write("# FONLL heavy-flavor meson grid\n")
        fh.write(f"# quark: {quark}\n# heavy_quark_mass_GeV: {mass:.8g}\n")
        fh.write(f"# scale: mu0=sqrt(m^2+pT^2); ffact(muF)={muF:.8g}, fren(muR)={muR:.8g}\n")
        fh.write("# pdf: NNPDF40_nlo_as_01180\n")
        fh.write(f"# lhapdf_id: {331700 + member}\n# lhapdf_member: {member}\n")
        fh.write(f"# variation_kind: {'central' if tag == 'central' else tag.split('_')[0]}\n")
        fh.write(f"# variation_tag: {tag}\n")
        fh.write(f"# trapezoid_integral_pb_y-3to3_pt0to50: {combine.trapz2(values):.12e}\n")
        fh.write("# pT y dsigma/dpT/dy\n")
        for pt, y, v in zip(ref_pt, ref_y, values):
            fh.write(f"{pt:.8g} {y:.8g} {v:.12e}\n")
    return {"quark": quark, "variation_kind": "central" if tag == "central" else tag.split("_")[0],
            "variation_tag": tag, "muR": muR, "muF": muF, "lhapdf_id": 331700 + member, "lhapdf_member": member,
            "heavy_quark_mass_GeV": mass, "path": str(path), "sha256": gen.sha256_file(path), "rows": 10000,
            "trapezoid_integral_pb": combine.trapz2(values)}


@pytest.fixture
def synthetic_campaign(tmp_path):
    out = tmp_path / "output"
    out.mkdir()
    entries = []
    spec = [("central", 0, 1.0, 1.0, 4.75, 1.0), ("scale_muR2p0_muF2p0", 0, 2.0, 2.0, 4.75, 0.9),
            ("scale_muR0p5_muF0p5", 0, 0.5, 0.5, 4.75, 1.2), ("pdf_0001", 1, 1.0, 1.0, 4.75, 1.05),
            ("pdf_0002", 2, 1.0, 1.0, 4.75, 0.95), ("mass_4p5", 0, 1.0, 1.0, 4.5, 1.1),
            ("mass_5p0", 0, 1.0, 1.0, 5.0, 0.92)]
    for tag, member, muR, muF, mass, scale in spec:
        entries.append(_write_grid(out / gen.grid_filename("nlo", tag, "bottom"), "bottom", tag, member, muR, muF,
                                   mass, scale))
    (out / "variation_manifest.json").write_text(json.dumps({"pdf_set": "NNPDF40_nlo_as_01180", "grids": entries}))
    return gen.Workspace(fonll=tmp_path / "fonll", out=out, run=tmp_path / "run", logs=tmp_path / "logs")


def test_partition_entries_splits_the_axes(synthetic_campaign):
    entries = json.loads((synthetic_campaign.out / "variation_manifest.json").read_text())["grids"]
    parts = combine.partition_entries(entries, "bottom")
    assert parts["central"]["variation_tag"] == "central"
    assert [e["variation_tag"] for e in parts["scale"]] == ["scale_muR0p5_muF0p5", "central", "scale_muR2p0_muF2p0"]
    assert [e["lhapdf_member"] for e in parts["pdf"]] == [0, 1, 2]
    assert [e["heavy_quark_mass_GeV"] for e in parts["mass"]] == [4.5, 4.75, 5.0]
    with pytest.raises(ValueError, match="no central grid"):
        combine.partition_entries(entries, "charm")


def test_combine_variations_writes_the_envelopes(synthetic_campaign):
    combine.combine_variations(synthetic_campaign, ["bottom", "charm"])
    env = synthetic_campaign.out / "envelopes"
    manifest = json.loads((env / "envelope_manifest.json").read_text())
    bands = [g["band"] for g in manifest["envelopes"]["bottom"]]
    assert bands == ["central", "scale_up", "scale_dn", "pdf_mean", "pdf_up", "pdf_dn", "mass_up", "mass_dn",
                     "combined_up", "combined_dn"]
    ref_pt, ref_y = combine.reference_grid_coords()
    central = combine.load_grid_column(Path(manifest["envelopes"]["bottom"][0]["path"]), ref_pt, ref_y)
    by_band = {g["band"]: combine.load_grid_column(Path(g["path"]), ref_pt, ref_y)
               for g in manifest["envelopes"]["bottom"]}
    assert np.allclose(by_band["scale_up"], 1.2 * central) and np.allclose(by_band["scale_dn"], 0.9 * central)
    assert np.allclose(by_band["mass_up"], 1.1 * central) and np.allclose(by_band["mass_dn"], 0.92 * central)
    sigma = np.std([1.05, 0.95], ddof=1)
    assert np.allclose(by_band["pdf_up"], central * (1 + sigma))
    quad = np.sqrt((0.2) ** 2 + sigma ** 2 + (0.1) ** 2)
    assert np.allclose(by_band["combined_up"], central * (1 + quad))
    # Every envelope carries the marker the production sampler refuses.
    assert "# envelope_band: combined_up" in Path(manifest["envelopes"]["bottom"][-2]["path"]).read_text()


def test_snapshot_rebuilds_an_equivalent_manifest(synthetic_campaign, monkeypatch):
    monkeypatch.setattr(snapshot, "STABLE_AGE_S", 0.0)
    original = json.loads((synthetic_campaign.out / "variation_manifest.json").read_text())["grids"]
    snapshot.snapshot(synthetic_campaign)
    rebuilt = json.loads((synthetic_campaign.out / "variation_manifest.json").read_text())
    assert rebuilt["note"].startswith("on-disk snapshot")
    assert rebuilt["pdf_members"] == [0, 1, 2]
    keys = ("quark", "variation_kind", "variation_tag", "muR", "muF", "lhapdf_id", "lhapdf_member",
            "heavy_quark_mass_GeV", "path", "sha256", "rows")
    assert [{k: e[k] for k in keys} for e in rebuilt["grids"]] == \
           [{k: e[k] for k in keys} for e in sorted(original, key=lambda e: (e["quark"], e["variation_tag"]))]


# ---------------------------------------------------------------- install --

def test_pins_name_every_shipped_diff():
    pins = json.loads((PATCHES / "PINS.json").read_text())
    diffs = {p.name for p in PATCHES.glob("*.diff")}
    assert {f["diff"] for f in pins["files"].values()} == diffs
    assert set(pins["files"]) == {"main/Makefile", "misc1/Makefile", "misc1/fonllgrid.f", "misc1/fragmfonll.f"}
    for rel, pin in pins["files"].items():
        text = (PATCHES / pin["diff"]).read_text()
        assert f"--- a/{rel}" in text and f"+++ b/{rel}" in text
        assert len(pin["sha256_pristine"]) == len(pin["sha256_patched"]) == 64


def test_bcfy_fragmentation_is_in_the_shipped_diff():
    text = (PATCHES / "misc1_fragmfonll.f.diff").read_text()
    assert "+      elseif(ifrag.eq.4) then" in text
    assert "+      elseif(ifrag.eq.5.or.ifrag.eq.8) then" in text
    assert text.count("+      parameter(nymx=120,nptmx=250,maxfiles=10)") == 3
    assert "+      parameter (nymx=120,nptmx=250)" in text


def test_install_status_classifies_files(tmp_path, monkeypatch):
    pristine = b"pristine\n"
    patched = b"patched\n"
    pins = {"fonll": {"version": "x", "mirror_commit": "abc"},
            "files": {"misc1/f.f": {"diff": "f.diff", "sha256_pristine": hashlib.sha256(pristine).hexdigest(),
                                    "sha256_patched": hashlib.sha256(patched).hexdigest()}}}
    monkeypatch.setattr(install, "PINS", pins)
    tree = tmp_path / "fonll"
    (tree / "misc1").mkdir(parents=True)
    assert install.status(tree) == {"misc1/f.f": "missing"}
    (tree / "misc1" / "f.f").write_bytes(pristine)
    assert install.status(tree) == {"misc1/f.f": "pristine"}
    (tree / "misc1" / "f.f").write_bytes(b"other\n")
    assert install.status(tree) == {"misc1/f.f": "unknown"}
    with pytest.raises(SystemExit, match="unknown"):
        install.apply_patches(tree)


def test_install_applies_a_diff_and_checks_the_result(tmp_path, monkeypatch):
    if subprocess.run(["patch", "--version"], capture_output=True).returncode != 0:
        pytest.skip("patch(1) not available")
    pristine = "line one\nline two\n"
    patched = "line one\nline 2\n"
    patch_dir = tmp_path / "patches"
    patch_dir.mkdir()
    (patch_dir / "f.diff").write_text("--- a/misc1/f.f\n+++ b/misc1/f.f\n@@ -1,2 +1,2 @@\n line one\n-line two\n+line 2\n")
    pins = {"fonll": {"version": "x", "mirror_commit": "abc"},
            "files": {"misc1/f.f": {"diff": "f.diff", "sha256_pristine": hashlib.sha256(pristine.encode()).hexdigest(),
                                    "sha256_patched": hashlib.sha256(patched.encode()).hexdigest()}}}
    monkeypatch.setattr(install, "PINS", pins)
    monkeypatch.setattr(install, "PATCHES", patch_dir)
    tree = tmp_path / "fonll"
    (tree / "misc1").mkdir(parents=True)
    (tree / "misc1" / "f.f").write_text(pristine)
    install.apply_patches(tree)
    assert (tree / "misc1" / "f.f").read_text() == patched
    assert install.status(tree) == {"misc1/f.f": "patched"}
    install.apply_patches(tree)      # idempotent
