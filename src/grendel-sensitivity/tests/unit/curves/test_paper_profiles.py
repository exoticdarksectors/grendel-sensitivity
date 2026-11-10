"""The manuscript profiles select layers from the diagnostic configuration
without mutating it."""
from pathlib import Path

from curves.plot import bc_panel
from curves.plot.paper import PAPER_PROJECTION_DRAWERS, paper_config


def _filenames(entries):
    return {Path(entry["path"]).name for entry in entries}


def test_bc4_paper_profile_keeps_only_direct_comparators():
    config = paper_config("bc4")
    assert _filenames(config["projections"]) == {"ship_bb_bc4.dat", "codexb_300fb.dat", "faser2_bc4.dat"}
    assert config["extra_excluded"] == []


def test_bc5_paper_profile_draws_only_the_projection_with_our_lifetime():
    config = paper_config("bc5")
    # SHiP uses our scalar lifetime; CODEX-b, MATHUSLA40 and ANUBIS assume shorter ones.
    assert _filenames(config["projections"]) == {"ship_pbc_bc5.dat"}
    # The PBC BC5 figure's existing exclusion is mixing-only and shared with BC4;
    # the ATLAS/CMS h -> SS bands of arXiv:2509.02564 are added to it.
    assert config["excluded"] == paper_config("bc4")["excluded"]
    assert _filenames(config["extra_excluded"]) == {"lhc_hss_hadronic_2402.01898.dat", "lhc_hss_mumu_2112.13769.dat",
                                                    "lhc_hss_mumu_2402.14491.dat", "lhc_hss_bb_llp_summaries.dat"}


def test_bc10_paper_profile_has_one_convention_and_compact_provenance():
    config = paper_config("bc10")
    # no CODEX-b: its published upper edge is out of reach of the B-decay production it assumes
    assert _filenames(config["projections"]) == {"ship_senscalc.dat", "faser2_senscalc_fpf.dat"}
    assert _filenames(config["extra_excluded"]) == {"beamdumps_past_alp2.dat", "kaon_charm_union.dat",
                                                    "kl_pi0ll.dat", "kpix_upper.dat", "na62_kpia_excluded.dat",
                                                    "e949_kpi_inv.dat", "babar.dat", "bs_mumu.dat", "lhcb_kmumu.dat", "lhcb_kstmumu.dat",
                                                    "na62_dump_2025.dat", "na62_kpix.dat", "na62_pi0_inv.dat"}
    assert config["secondary_ylabel"] is None


def test_bc10_bands_cover_the_unscanned_pole_windows():
    from grendel.models.alp.model import LIGHT_MESON_RESONANCE_WINDOWS
    table = Path(__file__).resolve().parents[3] / "curves" / "data" / "grendel_v2" / "bc10" / "sensitivity.csv"
    bands = bc_panel.unsupported_bands(table, bc_panel.BENCHMARKS["bc10"]["unsupported_labels"])
    assert sorted((lo, hi) for lo, hi, _ in bands) == [(0.52, 0.56), (0.94, 0.98)]   # computed neighbours
    for lo, hi, name in LIGHT_MESON_RESONANCE_WINDOWS:
        if name != "pi0":
            assert any(b_lo <= lo and hi <= b_hi for b_lo, b_hi, _ in bands)


def test_paper_profiles_do_not_mutate_the_diagnostic_configuration():
    paper_config("bc4")
    paper_config("bc5")
    paper_config("bc10")
    assert len(bc_panel.BENCHMARKS["bc4"]["projections"]) == 4
    assert len(bc_panel.BENCHMARKS["bc5"]["projections"]) == 4
    assert len(bc_panel.BENCHMARKS["bc10"]["extra_excluded"]) > 1
    assert bc_panel.BENCHMARKS["bc10"]["secondary_ylabel"]


def test_hnl_paper_profile_draws_ship_faser2_codexb_and_anubis():
    assert {drawer.__name__ for drawer in PAPER_PROJECTION_DRAWERS} == {
        "plot_ship_projection", "plot_faser2_projection", "plot_codexb_projection",
        "plot_anubis_projection"}


def test_every_reference_curve_the_profiles_name_is_tracked():
    for benchmark in ("bc4", "bc5", "bc10"):
        config = paper_config(benchmark)
        for entry in [*config["projections"], *config["extra_excluded"]]:
            assert Path(entry["path"]).is_file(), entry["path"]
        assert Path(config["excluded"]).is_file()


def test_hnlimits_data_files_resolve_regardless_of_letter_case(tmp_path):
    from curves.plot.hnl_panel import data_file
    (tmp_path / "ATLAS_2024").mkdir()
    (tmp_path / "ATLAS_2024" / "Ue4.dat").write_text("1 2\n")
    assert data_file("ATLAS_2024/Ue4.dat", tmp_path).read_text() == "1 2\n"
    assert data_file("atlas_2024/UE4.DAT", tmp_path).read_text() == "1 2\n"
    assert not data_file("atlas_2024/missing.dat", tmp_path).exists()


def _experiment(entry):
    from curves.plot.common import PROJECTION_STYLES
    label = (entry.get("label") or Path(entry["path"]).name).lower()
    return next((name for name in PROJECTION_STYLES if label.startswith(name.lower())), None)


def test_every_experiment_has_one_broken_style_on_every_panel():
    from curves.plot.common import PROJECTION_STYLES
    styles = list(PROJECTION_STYLES.values())
    assert len({s["color"] for s in styles}) == len(styles)
    assert len({s["ls"] for s in styles}) == len(styles)
    assert all(s["ls"] not in ("-", "solid") for s in styles)      # GRENDEL alone is solid
    for benchmark in ("bc4", "bc5", "bc10"):
        for entry in paper_config(benchmark)["projections"]:
            name = _experiment(entry)
            assert name is not None, entry["path"]
            assert (entry["color"], entry["ls"]) == (PROJECTION_STYLES[name]["color"], PROJECTION_STYLES[name]["ls"])


def test_grendel_threshold_contour_is_solid():
    import inspect
    from curves.plot import hnl_panel
    for function in (bc_panel.plot_grendel_threshold, hnl_panel.plot_grendel_threshold):
        assert inspect.signature(function).parameters["ls"].default == "-"


def test_paper_profile_frames_and_labels():
    bc4, bc5 = paper_config("bc4"), paper_config("bc5")
    assert bc4["x_range"][1] == 4.0
    assert all(text is None for _, text in bc5["x_markers"])          # m_h/2: the line, no label
    for name in ("bc4", "bc5", "bc10"):
        ship = next(e for e in paper_config(name)["projections"] if e["path"].name.startswith("ship_"))
        assert ship["label"].startswith("SHiP ($6\\times10^{20}$ PoT)")
