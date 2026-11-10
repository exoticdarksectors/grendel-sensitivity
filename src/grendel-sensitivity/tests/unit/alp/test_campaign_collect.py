"""BC10 collect refuses variations produced from different code states before
it writes any output."""
import pandas as pd
import pytest

from grendel.band.alp import campaign


def test_mixed_code_states_are_refused_before_any_output(tmp_path, monkeypatch):
    registry = [{"code": {"commit": "a"}}, {"code": {"commit": "b"}}]
    monkeypatch.setattr(campaign, "load_completed", lambda scratch, variations: (pd.DataFrame(), registry))
    monkeypatch.setattr(campaign, "combine_band", lambda raw, central: pytest.fail("combined before the check"))
    out = tmp_path / "band"
    out.mkdir()
    with pytest.raises(ValueError, match="different code states"):
        campaign.collect(tmp_path / "scratch", [], tmp_path / "central.csv", out)
    assert list(out.iterdir()) == []
