"""The ATLAS/CMS h -> S S bands of the BC5 panel (arXiv:2509.02564)."""
from pathlib import Path

import numpy as np

from curves.plot import paper
from curves.plot.common import load_polygon_segments
from curves.tools import lhc_bc5


def test_the_four_searches_are_bands_on_the_paper_panel():
    config = paper.paper_config("bc5")
    files = [Path(e["path"]) for e in config["extra_excluded"]]
    assert {f.name for f in files} == {f"{stem}.dat" for stem, _ in lhc_bc5.SEARCHES.values()}
    assert all(e["mode"] == "band" for e in config["extra_excluded"])
    for f in files:
        chains = load_polygon_segments(f)
        assert len(chains) % 2 == 0
        for lower, upper in zip(chains[0::2], chains[1::2]):
            assert np.array_equal(lower[:, 0], upper[:, 0]) and np.all(np.diff(lower[:, 0]) > 0)
            assert np.all((lower[:, 1] > 0) & (lower[:, 1] < upper[:, 1]))


def test_the_dimuon_band_is_split_at_the_vetoed_resonances():
    rows = [(m, 1e-10, 1e-8) for m in (2, 2.91, 3.89, 5, 8.99, 10.78, 12)]
    split = lhc_bc5.pieces(rows, lhc_bc5.GAPS["arXiv:2112.13769"])
    assert [[r[0] for r in piece] for piece in split] == [[2, 2.91], [3.89, 5, 8.99], [10.78, 12]]


def test_only_the_one_percent_scenario_is_read():
    text = ("# Scenario: $\\mathcal{B}(h\\rightarrow SS) = 1\\%$\n"
            "# Label: $\\mathrm{S}\\to \\mu\\mu$ (arXiv:2402.14491)\n"
            "m(S) [GeV], min_excluded, max_excluded\n10, 1e-10, 1e-8\n"
            "# Scenario: $\\mathcal{B}(h\\rightarrow SS) = 0.1\\%$\n"
            "# Label: $\\mathrm{S}\\to \\mu\\mu$ (arXiv:2402.14491)\n"
            "m(S) [GeV], min_excluded, max_excluded\n10, 2e-10, 1e-9\n")
    assert list(lhc_bc5.parse(text).values()) == [[(10.0, 1e-10, 1e-8)]]
