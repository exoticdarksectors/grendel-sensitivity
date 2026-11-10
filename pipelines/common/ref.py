"""Where the validation steps find their reference material (the Python side of pipelines/common/ref.sh).

No production step reads any of it. Three kinds:
    PUBLISHED      the published tables: src/grendel-sensitivity/curves/data (in this repository)
    PAPER_FIGURES  the figure files of the paper: reference/paper_figures (in this repository)
    REF            the outputs of the original runs (campaign records, event samples, templates, ...), $GRENDEL_REF:
                   optional and not in this repository (reference/README.md; GRENDEL_REF=none: without it).
                   Without it every comparison that needs it prints a SKIP line, collected in
                   work/skipped_validations.txt.

    sys.path.insert(0, str(Path(os.environ["REPRO"]) / "pipelines/common")); from ref import REF, have_ref
    if have_ref("BC10 templates vs the stored bundles", REF / "grendel/campaigns/..."): ...
    skip_unless("...", path)    # in a script that is nothing but a validation: exits 0 after the SKIP line
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

REPRO = Path(os.environ.get("REPRO") or Path(__file__).resolve().parents[2])
PUBLISHED = Path(os.environ.get("GRENDEL_PUBLISHED") or REPRO / "src/grendel-sensitivity/curves/data")
PAPER_FIGURES = Path(os.environ.get("GRENDEL_PAPER_FIGURES") or REPRO / "reference/paper_figures")
_root = os.environ.get("GRENDEL_REF", "none")
REF = Path(_root) if _root not in ("", "none") else Path("/nonexistent/GRENDEL_REF")  # "REF / ..." never exists without a tree


def have_ref(label: str, *paths: Path | str) -> bool:
    """True if every path exists; otherwise print 'SKIP label: ...', log it and return False."""
    missing = [Path(p) for p in paths if not Path(p).exists()]
    if not missing:
        return True
    first = str(missing[0])
    if first.startswith(str(REF) + "/"):
        first = first[len(str(REF)) + 1:]
    what = first + (f" and {len(missing) - 1} more" if len(missing) > 1 else "")
    msg = f"SKIP {label}: reference data not available ({what}; see reference/README.md)"
    print(msg, file=sys.stderr, flush=True)
    log = REPRO / "work/skipped_validations.txt"
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("a") as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\t{msg}\n")
    return False


def skip_unless(label: str, *paths: Path | str) -> None:
    """For a script that only validates: exit 0 (after the SKIP line) when its reference data is missing."""
    if not have_ref(label, *paths):
        sys.exit(0)
