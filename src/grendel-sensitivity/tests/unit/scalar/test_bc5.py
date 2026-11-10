"""BC5: the quartic coupling, the high-mass decay model, the coupling-power
rows of the scan and the production channels."""
import gzip

import numpy as np
import pytest

from grendel.io.vectors import load_combined_csv, write_llp_csv
from grendel.models.scalar import model, production_bc5
from grendel.models.scalar.higgs_pool import M_HIGGS, sample_higgs_4vectors
from grendel.reco.acceptance import scan_u2, signal_contribution_diagnostics


# --- the quartic coupling ---------------------------------------------------

def test_br_h_to_SS_matches_boiarska_eq_2_15():
    # alpha = 1 GeV gives BR(h -> SS) ~ 2.0e-2 (arXiv:1904.10447 eq. 2.15)
    assert model.br_h_to_SS(1.0, alpha=1.0) == pytest.approx(0.020, rel=0.03)
    alpha = model.alpha_quartic(1.0, br_hss=0.01)
    assert model.br_h_to_SS(1.0, alpha) == pytest.approx(0.01, rel=1e-9)
    assert model.alpha_quartic(70.0) == 0.0          # closed above m_h / 2


def test_B_to_K_SS_reproduces_table_4_and_senscalc():
    # Boiarska Table 4: BR(B+ -> K+ S S) = 1.4e-10 at m_S = 0, alpha = 1 GeV;
    # SensCalc's table starts at 1.357e-10.
    m, br_table = model._senscalc_table("BplustoKplusSS")
    assert model.br_B_to_K_SS(1e-3, alpha=1.0) == pytest.approx(br_table[0], rel=0.03)
    assert model.br_B_to_K_SS(1.0, alpha=1.0) == pytest.approx(
        np.interp(1.0, m, br_table), rel=0.05)
    assert model.br_B_to_K_SS(2.0, alpha=1.0) < model.br_B_to_K_SS(1.0, alpha=1.0)
    assert model.br_B_to_K_SS(2.5, alpha=1.0) == 0.0  # 2 m_S > m_B - m_K


def test_quartic_rates_scale_as_alpha_squared_and_close():
    for fn in (model.br_b_to_Xs_SS, model.br_Bs_to_SS):
        assert fn(1.0, alpha=0.5) == pytest.approx(0.25 * fn(1.0, alpha=1.0))
        assert fn(1.0, alpha=1.0) > 0.0
        assert fn(6.0, alpha=1.0) == 0.0
    # inclusive > exclusive K channel; Lambda_b enters through its lifetime
    assert model.br_b_to_Xs_SS(0.5) > model.br_B_to_K_SS(0.5)
    ratio = model.br_b_to_Xs_SS(0.5, parent="Lambda_b") / model.br_b_to_Xs_SS(0.5, parent="B+")
    assert ratio == pytest.approx(model.TAU_LAMBDA_B / model.TAU_BPLUS)


def test_q2_sampler_stays_inside_the_dalitz_range():
    rng = np.random.default_rng(3)
    m_S, m_B, m_K = 0.7, model.M_BPLUS, model.M_KPLUS
    q2 = model.sample_q2_B_to_K_SS(5000, m_S, rng, m_B, m_K)
    assert len(q2) == 5000
    assert q2.min() >= 4 * m_S ** 2 and q2.max() <= (m_B - m_K) ** 2
    assert len(model.sample_q2_B_to_K_SS(10, 3.0, rng, m_B, m_K)) == 0


# --- the decay model above the Winkler digitisation ---------------------------

def test_widths_continue_across_the_winkler_end():
    below, above = model.ctau_sin2theta1(7.45), model.ctau_sin2theta1(7.55)
    assert 0.6 < above / below < 1.5
    for m_S in (10.0, 20.0, 40.0, 62.0):
        br = model.branching_ratios(m_S)
        assert sum(br.values()) == pytest.approx(1.0)
        assert set(br) == set(model.branching_ratios(3.0))
    assert model.branching_ratios(20.0)["bb"] > 0.7
    assert model.branching_ratios(7.0)["bb"] == 0.0
    assert model.ctau_sin2theta1(20.0) == pytest.approx(4.33e-13, rel=0.05)
    spectator = model.partial_widths(20.0, scheme="chpt_spectator")
    assert spectator["bb"] > spectator["cc"] > 0.0


# --- coupling-power rows in the scan ------------------------------------------

def _toy_mc(rng, n_ev=40, n_samples=8):
    d = rng.uniform(0.0, 5.0, (n_ev, n_samples))
    passed = rng.random((n_ev, n_samples)) < 0.6
    path = np.full(n_ev, 5.0)
    weight = rng.uniform(0.5, 2.0, n_ev)
    bg = rng.uniform(2.0, 20.0, n_ev)
    return d, passed, path, weight, bg


def test_scan_u2_coupling_power():
    d, passed, path, weight, bg = _toy_mc(np.random.default_rng(5))
    grid = np.logspace(-10, -4, 7)
    _, plain = scan_u2(d, passed, path, weight, bg, 1e-6, 3e6, grid)
    _, ones = scan_u2(d, passed, path, weight, bg, 1e-6, 3e6, grid, coupling_power=np.ones(len(weight)))
    _, zeros = scan_u2(d, passed, path, weight, bg, 1e-6, 3e6, grid, coupling_power=np.zeros(len(weight)))
    np.testing.assert_allclose(plain, ones, rtol=1e-12)   # same sum, u2 applied per row
    np.testing.assert_allclose(zeros, plain / grid, rtol=1e-12)
    power = (np.arange(len(weight)) % 2).astype(float)
    _, mixed = scan_u2(d, passed, path, weight, bg, 1e-6, 3e6, grid, coupling_power=power)
    _, part1 = scan_u2(d, passed, path, weight * power, bg, 1e-6, 3e6, grid)
    _, part0 = scan_u2(d, passed, path, weight * (1 - power), bg, 1e-6, 3e6, grid)
    np.testing.assert_allclose(mixed, part1 + part0 / grid, rtol=1e-12)
    diag = signal_contribution_diagnostics(d, passed, path, weight, bg, 1e-6, 1e-7, coupling_power=power)
    assert 0 < diag["event_ess"] <= len(weight)


def test_load_combined_csv_reads_the_optional_sixth_column(tmp_path):
    w, E, px, py, pz = np.ones(3), np.full(3, 5.0), np.zeros(3), np.zeros(3), np.full(3, 4.0)
    write_llp_csv(tmp_path / "five.csv", w, E, px, py, pz)
    write_llp_csv(tmp_path / "six.csv", w, E, px, py, pz, coupling_power=[1, 0, 0])
    five = load_combined_csv(tmp_path / "five.csv", 3.0)
    six = load_combined_csv(tmp_path / "six.csv", 3.0)
    np.testing.assert_array_equal(five["coupling_power"], [1.0, 1.0, 1.0])
    np.testing.assert_array_equal(six["coupling_power"], [1.0, 0.0, 0.0])
    np.testing.assert_array_equal(five["E"], six["E"])


# --- production -----------------------------------------------------------------

@pytest.fixture
def toy_higgs_pool(tmp_path):
    path = tmp_path / "pool.csv.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("# toy\n")
        fh.write("ggf,0,20,-1,1,300\nggf,20,60,-3,-1,100\nvbf,40,80,1,3,50\n")
    return path


def test_higgs_sampler_follows_the_histogram(toy_higgs_pool):
    rng = np.random.default_rng(1)
    higgs = sample_higgs_4vectors(20000, rng, sigma_pb={"ggf": 3.0, "vbf": 1.0}, path=toy_higgs_pool)
    assert higgs["sigma_pb"] == pytest.approx(4.0)
    assert np.mean(higgs["mode"] == "vbf") == pytest.approx(0.25, abs=0.02)
    in_first = (higgs["pt"] < 20) & (np.abs(higgs["y"]) < 1)
    assert in_first.mean() == pytest.approx(0.75 * 0.75, abs=0.02)
    m2 = higgs["E"] ** 2 - higgs["px"] ** 2 - higgs["py"] ** 2 - higgs["pz"] ** 2
    np.testing.assert_allclose(np.sqrt(m2), M_HIGGS, rtol=1e-9)


def test_bc5_rows_carry_the_right_coupling_power(toy_higgs_pool):
    rngs = production_bc5.make_streams(7)
    from grendel.production.fonll.meson_sampler import sample_meson_4vectors
    pool = sample_meson_4vectors(400, "bottom", rng=rngs[0])
    higgs = sample_higgs_4vectors(300, rngs[2], sigma_pb={"ggf": 50.0}, path=toy_higgs_pool)
    (w, E, px, py, pz, power, channel), counts = production_bc5.generate_bc5_4vectors(
        1.0, rngs, pool=pool, higgs=higgs, sigma_bottom=5e8, br_hss=0.01)
    assert counts["hSS"] == 600 and counts["BSS"] > 0 and counts["mixing"] > 0
    assert len(w) == sum(counts.values()) == len(power)
    assert set(power[:counts["mixing"]]) == {1.0} and set(power[counts["mixing"]:]) == {0.0}
    assert list(np.unique(channel, return_counts=True)[1]) == [counts[name] for name in production_bc5.CHANNELS]
    hss = slice(counts["mixing"], counts["mixing"] + counts["hSS"])
    assert w[hss].sum() == pytest.approx(2 * 50.0 * 0.01)
    mass2 = E ** 2 - px ** 2 - py ** 2 - pz ** 2
    np.testing.assert_allclose(np.sqrt(mass2), 1.0, rtol=1e-6)
    # above the B ceiling only h -> S S remains; above m_h/2 nothing
    (_, _, _, _, _, power20, _), counts20 = production_bc5.generate_bc5_4vectors(
        20.0, rngs, pool=pool, higgs=higgs, sigma_bottom=5e8)
    assert counts20 == {"mixing": 0, "hSS": 600, "BSS": 0}
    _, counts70 = production_bc5.generate_bc5_4vectors(70.0, rngs, pool=pool, higgs=higgs, sigma_bottom=5e8)
    assert sum(counts70.values()) == 0
    # BR(h -> SS) = 0 and no quartic B: exactly the BC4 rows
    (w0, *_, p0, _), counts0 = production_bc5.generate_bc5_4vectors(
        1.0, production_bc5.make_streams(7), pool=pool, higgs=None, sigma_bottom=5e8, br_hss=0.0, quartic_b=False)
    assert counts0["hSS"] == 0 and counts0["BSS"] == 0 and set(p0) == {1.0}
