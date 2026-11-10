"""BC4 --templates-dir defaults to the analytic decay engine, not to an
environment variable (the option is shared with the HNL and BC10 scans)."""
import pytest

from grendel.models.hnl import scan as hnl_scan
from grendel.models.scalar import scan as scalar_scan


def _help(main, capsys):
    with pytest.raises(SystemExit):
        main(["--help"])
    return " ".join(capsys.readouterr().out.split())


def test_templates_dir_help_names_each_default(capsys):
    bc4 = _help(scalar_scan.main, capsys)
    assert "default: none, the analytic decay engine" in bc4 and "GRENDEL_BC4_TEMPLATES_DIR" not in bc4
    assert "decay templates (default: $GRENDEL_HNL_TEMPLATES_DIR)" in _help(hnl_scan.main, capsys)
