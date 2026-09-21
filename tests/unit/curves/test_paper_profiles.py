"""The manuscript profiles select layers from the diagnostic configuration
without mutating it."""
from pathlib import Path

from curves.plot import bc_panel
from curves.plot.paper import PAPER_PROJECTION_DRAWERS, paper_config


def _filenames(entries):
    return {Path(entry["path"]).name for entry in entries}


def test_bc4_paper_profile_keeps_only_direct_comparators():
    config = paper_config("bc4")
    assert _filenames(config["projections"]) == {"ship_bc4.dat", "codexb_300fb.dat"}
    assert config["extra_excluded"] == []


def test_bc10_paper_profile_has_one_convention_and_compact_provenance():
    config = paper_config("bc10")
    assert _filenames(config["projections"]) == {"ship_alp2.dat"}
    assert _filenames(config["extra_excluded"]) == {"beamdumps_past_alp2.dat"}
    assert config["secondary_ylabel"] is None


def test_paper_profiles_do_not_mutate_the_diagnostic_configuration():
    paper_config("bc4")
    paper_config("bc10")
    assert len(bc_panel.BENCHMARKS["bc4"]["projections"]) == 3
    assert len(bc_panel.BENCHMARKS["bc10"]["extra_excluded"]) > 1
    assert bc_panel.BENCHMARKS["bc10"]["secondary_ylabel"]


def test_hnl_paper_profile_draws_ship_faser2_and_codexb():
    assert {drawer.__name__ for drawer in PAPER_PROJECTION_DRAWERS} == {
        "plot_ship_projection", "plot_faser2_projection", "plot_codexb_projection"}


def test_every_reference_curve_the_profiles_name_is_tracked():
    for benchmark in ("bc4", "bc10"):
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
