"""The BC4 campaign config records what a run produces: the per-mass scalar
event count (Lambda_b closes at m_Lambda_b - m_Lambda = 4.50 GeV) and every
source file that determines the result."""
import os
import subprocess
import sys
from pathlib import Path

from grendel.band.scalar import campaign
from grendel.io.paths import repo_root


def test_expected_events_follow_the_species_ceilings(tmp_path, monkeypatch):
    monkeypatch.setattr(campaign, "ray_backend_provenance", lambda: {"implementation": "stub"})
    grid = tmp_path / "central.dat"
    grid.write_text("grid\n")
    variation = {"name": "central", "axis": "central", "width_scheme": "winkler", "path": grid, "sha256": "0" * 64}
    config = campaign.campaign_config(variation, [0.5, 4.5, 4.6, 4.7], 10, 3, 42)
    assert config["n_scalar_events_expected"] == [40, 40, 30, 30]


def test_code_inputs_cover_the_result_determining_files():
    for rel in ("grendel/models/scalar/data/winkler_widths.csv", "grendel/production/constants.py",
                "grendel/production/fonll/fonll_parser.py", "grendel/production/fonll/meson_sampler.py",
                "grendel/production/decay_engine/kinematics.py", "grendel/constants.py",
                "grendel/models/hnl/spec.py", "grendel/models/scalar/templates_exhad.py"):
        assert rel in campaign.CODE_INPUTS, rel
    assert all((repo_root() / rel).is_file() for rel in campaign.CODE_INPUTS)


def test_code_inputs_cover_every_module_the_campaign_imports():
    # A fresh interpreter, so only the campaign's own imports are counted. The
    # band bookkeeping and the path/atomic-write helpers do not change a result.
    probe = ("import sys, grendel.band.scalar.campaign, grendel.models.scalar.spec, grendel.scan\n"
             "for n, m in list(sys.modules.items()):\n"
             "    if n.startswith('grendel') and getattr(m, '__file__', None):\n"
             "        print('MODULE', m.__file__)")
    out = subprocess.run([sys.executable, "-c", probe], cwd=repo_root(), capture_output=True, text=True,
                         check=True, env={**os.environ, "PYTHONPATH": str(repo_root())}).stdout
    imported = {Path(line[7:]).resolve().relative_to(repo_root()).as_posix()
                for line in out.splitlines() if line.startswith("MODULE ")}
    assert "grendel/scan.py" in imported
    exempt = {"grendel/io/atomic.py", "grendel/io/paths.py"}
    missing = {f for f in imported - exempt - set(campaign.CODE_INPUTS)
               if not f.endswith("__init__.py") and not f.startswith("grendel/band/")}
    assert not missing
