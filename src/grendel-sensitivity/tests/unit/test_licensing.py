"""The FairShip-derived files keep their LGPL notice.

Two files are adaptations of FairShip sources rather than callers of them, so
they carry FairShip's licence while the rest of the repository is
BSD-3-Clause. This pins that boundary: the notice cannot be dropped by an
edit, and no other ``.py`` or ``.conf`` file can acquire one without being documented.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

LGPL_FILES = (
    "grendel/models/hnl/templates/fairship_decay.py",
    "grendel/models/hnl/data/fairship_decay_selection.conf",
)
TAG = "SPDX-License-Identifier: LGPL-3.0-or-later"


def test_licence_text_is_shipped():
    text = (ROOT / "LICENSES" / "LGPL-3.0-or-later.txt").read_text()
    assert "GNU LESSER GENERAL PUBLIC LICENSE" in text
    assert "GNU GENERAL PUBLIC LICENSE" in text  # LGPLv3 incorporates GPLv3


def test_derived_files_carry_the_notice():
    for rel in LGPL_FILES:
        head = (ROOT / rel).read_text()[:400]
        assert TAG in head, rel
        assert "SHiP Collaboration" in head, rel


def test_the_notice_and_the_documentation_agree():
    documented = (ROOT / "THIRD_PARTY.md").read_text()
    for rel in LGPL_FILES:
        assert rel in documented, rel

    tagged = {
        str(path.relative_to(ROOT))
        for path in ROOT.rglob("*")
        if path.is_file()
        and path.suffix in {".py", ".conf"}
        and ".git" not in path.parts
        and "third_party" not in path.relative_to(ROOT).parts[:1]
        and path != Path(__file__)
        and TAG in path.read_text(errors="ignore")[:400]
    }
    assert tagged == set(LGPL_FILES)
