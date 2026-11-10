"""The comparison layers corrected after v1: BBN, CODEX-b, SHiP and FASER2 BC4, beam dumps."""
import hashlib
import math
import re
from pathlib import Path

import numpy as np
import pytest

from curves.plot import bc_panel, hnl_panel
from curves.plot.common import load_polygon_segments
from curves.plot.data import BC4_DATA, DATA, HNL_DATA
from curves.plot.paper import REGION_LABELS, paper_config
from curves.plot.results import GrendelResults

FLAVORS = ("Ue", "Umu", "Utau")


# ------------------------------------------------------------------- BBN --

@pytest.mark.parametrize("flavor", FLAVORS)
def test_bbn_layer_is_the_single_majorana_region_as_published(flavor):
    path = HNL_DATA / hnl_panel.COSMO_PER_FLAVOUR_FILES[flavor]
    assert "arXiv:2410.07343" in path.read_text()
    region = hnl_panel.load_region_file(path)
    record = {"id": "cosmo", "file_top": "Sabti_BBN/Ue4_top_Cosmo.dat",
              "file_bottom": "Sabti_BBN/Ue4_bottom_Cosmo.dat", "hnl_type": "Dirac", "CL": "2sigma", "type": "C"}
    scenario = next(s for s, c in hnl_panel.SCENARIOS.items() if c["flavor"] == flavor)
    (item,) = hnl_panel.processed_curves(scenario, hnl_panel.SCENARIOS[scenario], [record])
    assert item["category"] == "cosmology"
    assert np.array_equal(item["top"], region["top"]) and np.array_equal(item["bottom"], region["bottom"])
    top = region["top"]
    # the upper edge leaves the 1e-10 panel floor before the calculation stops
    assert top[-1, 1] < hnl_panel.Y_RANGE[0] < top[top[:, 0] >= 0.3][0, 1]


def test_bbn_legend_names_the_model():
    assert "single Majorana" in REGION_LABELS["cosmology"]


# --------------------------------------------------------------- CODEX-b --

def _codexb(flavor):
    rows = [line.split() for line in (HNL_DATA / "original_sources" / "CODEX-b" / f"CODEX-b_2019_{flavor}.dat")
            .read_text().splitlines() if line.strip() and not line.startswith("#")]
    return np.asarray(rows, dtype=float)


@pytest.mark.parametrize("flavor", FLAVORS)
def test_codexb_majorana_curve_is_open_at_the_frame(flavor):
    rows = _codexb(flavor)
    majorana, dirac = rows[np.isfinite(rows[:, 1])], rows[np.isfinite(rows[:, 2])]
    assert dirac[:, 2].max() == pytest.approx(1.0e-2)       # the EoI contour leaves its frame
    assert majorana[:, 1].max() < 2.0e-3                     # the Majorana curve has no frame vertex
    paired = rows[np.isfinite(rows[:, 1]) & np.isfinite(rows[:, 2])]
    ratio = paired[:, 1] / paired[:, 2]
    assert ratio.max() <= 1.0 and ratio.min() > 0.3
    assert ratio[0] == pytest.approx(1.0 / math.sqrt(2.0))  # lowest mass: lower edge, frame-limited column


def test_codexb_conversion_rescales_each_edge():
    fitz = pytest.importorskip("fitz")  # noqa: F841 - the tool module imports PyMuPDF
    from curves.tools.vector_codexb_hnl import majorana_rows
    # lower edge at m = 1, 2, 3; upper edge at 3, 2; frame vertex
    rows = [(1.0, 1e-7), (2.0, 1e-8), (3.0, 1e-7), (3.0, 3e-7), (2.0, 1e-4), (1.5, 1e-2)]
    out = majorana_rows(rows)
    assert [r[0] for r in out] == [r[0] for r in rows]                  # no vertex added or moved in mass
    assert all(out[i][1] == pytest.approx(rows[i][1] / math.sqrt(2.0)) for i in range(3))   # lower edge
    assert all(out[i][1] == pytest.approx(rows[i][1] / 2.0) for i in (3, 4))                # upper edge
    assert math.isnan(out[5][1]) and out[5][2] == 1e-2                  # frame vertex: Dirac only
    with pytest.raises(RuntimeError, match="cross"):                    # edges closer than sqrt2
        majorana_rows([(1.0, 1e-7), (2.0, 1e-7), (2.0, 1.2e-7)])


# ------------------------------------------------------------- SHiP BC4 --

def test_bc4_ship_and_faser2_are_the_briefing_book_curves_and_older_ship_traces_are_not_drawn():
    for name in ("ship_bb_bc4.dat", "faser2_bc4.dat"):
        header = (BC4_DATA / "projections" / name).read_text()
        assert "arXiv:2511.03883" in header and "Fig. 8.17" in header and "Br(h -> SS) = 0" in header
    assert "arXiv:2310.17726" in (BC4_DATA / "projections" / "ship_bc4.dat").read_text()
    assert (BC4_DATA / "projections" / "ship_bc4_2504_reference.dat").is_file()
    drawn = {Path(e["path"]).name for e in paper_config("bc4")["projections"]}
    assert drawn == {"ship_bb_bc4.dat", "codexb_300fb.dat", "faser2_bc4.dat"}


# ----------------------------------------------------------- beam dumps --

def test_bc10_beam_dumps_are_clipped_to_their_source_frame():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    config = paper_config("bc10")
    entries = {Path(e["path"]).name: e for e in config["extra_excluded"]}
    entry = entries["beamdumps_past_alp2.dat"]
    assert entry["x_clip"] == (0.6, 2.0)
    assert min(seg[:, 0].min() for seg in load_polygon_segments(entry["path"])) < 0.6   # data beneath the clip
    # the FIPs 2022 CHARM/kaon region hands over to it at 0.6 GeV
    union = entries["kaon_charm_union.dat"]
    assert union["x_clip"] == (0.01, 0.6)
    assert max(seg[:, 0].max() for seg in load_polygon_segments(union["path"])) > 0.6
    fig, ax = plt.subplots()
    bc_panel.plot_extra_excluded(ax, {**config, "extra_excluded": [entry, union]})
    for patch in ax.patches:
        assert patch.get_clip_path() is not None and patch.get_clip_box() is not None   # axes clip kept
    plt.close(fig)


# ------------------------------------------------------------ diagnostic --

def test_raw_hnlimits_copies_are_read_in_their_own_unit(tmp_path):
    from curves.plot.diagnostic import load_raw_record
    path = tmp_path / "x_Ue.dat"
    path.write_text("# mass_MeV  U2\n# bottom\n100  1e-6\n300  1e-8\n")
    bottom, top = load_raw_record(path)
    assert top is None and np.allclose(bottom[:, 0], [0.1, 0.3])


def test_diagnostic_regions_are_drawn_as_hnlimits_draws_them():
    from curves.plot.diagnostic import excluded_mask, region_mask
    m, u2 = np.array([0.5, 1.0, 2.0]), np.array([1e-9, 1e-6, 1e-3])
    bottom, top = np.array([[0.2, 1e-7], [3.0, 1e-7]]), np.array([[0.2, 1e-4], [3.0, 1e-4]])
    assert region_mask(bottom, None, m, u2)[:, 1].tolist() == [False, True, True]
    assert region_mask(bottom, top, m, u2)[:, 1].tolist() == [False, True, False]
    box = np.array([[0.8, 1e-7], [1.5, 1e-7], [1.5, 1e-5], [0.8, 1e-5]])
    assert region_mask(box, box, m, u2)[:, 1].tolist() == [False, True, False]
    # cosmology is the BBN region: excluded below its upper edge, not above it
    grid = np.array([0.5])
    top_05 = hnl_panel.log_interp(*hnl_panel.load_region_file(HNL_DATA / hnl_panel.COSMO_PER_FLAVOUR_FILES["Ue"])
                                  ["top"].T, grid)[0]
    mask = excluded_mask("Ue", grid, np.array([top_05 / 3, top_05 * 3]), category="cosmology")
    assert mask[:, 0].tolist() == [True, False]


# --------------------------------------------------- sources, v1 tables --

def test_hand_maintained_source_notes_are_left_alone():
    from curves.tools.original_sources import SOURCES
    flags = {s["exp"]: bool(s.get("hand_maintained")) for s in SOURCES}
    assert flags == {"SHiP": False, "MATHUSLA": False, "CODEX-b": True, "ANUBIS": True}


V1_MODELS = ("hnl", "bc4", "bc10")     # the benchmarks of arXiv:2609.00152v1
V2_MODELS = ("hnl", "bc4", "bc5", "bc10")


@pytest.mark.parametrize("release, models", [("grendel_v1", V1_MODELS), ("grendel_v2", V2_MODELS)])
def test_the_release_tables_are_tracked_with_their_hashes(release, models):
    readme = (DATA / release / "README.md").read_text()
    listed = dict(re.findall(r"\| `(\w+/sensitivity\w*\.csv)` \|.*`([0-9a-f]{64})` \|", readme))
    results = GrendelResults(DATA / release)
    assert len(listed) == 2 * len(models)
    for model in models:
        for path in (results.sensitivity(model), results.nsig10(model)):
            rel = f"{model}/{path.name}"
            assert hashlib.sha256(path.read_bytes()).hexdigest() == listed[rel], rel
