"""Read the stored (published) HNL four-vectors of the reference campaigns and split a stored combined file into its
channel blocks.

The published coarse-grid rows were scanned on asbuilt_20261014/hnl/llp_4vectors/<fl>/combined/mN_<label>.csv, the
concatenation (grendel.models.hnl.production.combine_channels.CHANNELS order)
    Bmeson, Dmeson, Bc, Bbaryon, tau, induced_tau, Kmeson, WZ
of the channel files. Dmeson, Bc and Kmeson are readable next to it; the Bmeson, Bbaryon, tau, induced_tau and WZ
files of the June v1 production are links to an HFS+-compressed tree that reads as empty on Linux, so they are
recovered from the combined file: every block length except the MadGraph prompt-tau one is fixed by the recipe
(n_pool, species split, 100k W/Z events), and the tau block is what remains.

The campaigns are part of the optional reference tree $GRENDEL_REF (pipelines/common/ref.py, reference/README.md);
without it the paths below do not exist and the validations that need them are skipped.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(os.environ["REPRO"]) / "pipelines/common"))
from ref import REF as _REF  # noqa: E402

REF = _REF / "grendel/campaigns"
ASBUILT_VECTORS = REF / "asbuilt_20261014/hnl/llp_4vectors"
FINE = REF / "hnl/hnl_fine_condor_20261014"
ORDER = ["Bmeson", "Dmeson", "Bc", "Bbaryon", "tau", "induced_tau", "Kmeson", "WZ"]


def label(m: float) -> str:
    return f"{m:.3f}".replace(".", "p")


def read_rows(path: Path) -> list[str]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with open(path) as f:
        return f.read().splitlines()


def split_stored_combined(flavor: str, mass: float, n_bmeson: int, n_bbaryon: int, n_wz: int = 100_000
                          ) -> dict[str, list[str]]:
    """Split the stored combined file into channel blocks (text rows, as written).

    n_bmeson, n_bbaryon, n_wz: block lengths fixed by the recipe (taken from a production with the same n_pool;
    the June induced-tau pool had 99 999 rows where the current code builds 100 000, so its length is not used).
    Dmeson, Bc and Kmeson lengths come from the readable stored channel files, which are checked against their
    blocks. What lies between Bbaryon and Kmeson is the prompt-tau block (one run of equal weights: every pool
    row carries sigma*K/N*BR) followed by the induced-tau block (one equal-weight run per tau source)."""
    rows = read_rows(ASBUILT_VECTORS / flavor / "combined" / f"mN_{label(mass)}.csv")
    own = {ch: read_rows(ASBUILT_VECTORS / flavor / ch / f"mN_{label(mass)}.csv") for ch in ("Dmeson", "Bc", "Kmeson")}
    blocks, i = {}, 0
    for ch, n in (("Bmeson", n_bmeson), ("Dmeson", len(own["Dmeson"])), ("Bc", len(own["Bc"])),
                  ("Bbaryon", n_bbaryon)):
        blocks[ch] = rows[i:i + n]
        i += n
    j = len(rows) - n_wz
    blocks["WZ"] = rows[j:]
    k = j - len(own["Kmeson"])
    blocks["Kmeson"] = rows[k:j]
    middle = rows[i:k]
    n_tau = 0
    if middle:
        w0 = middle[0].split(",", 1)[0]
        while n_tau < len(middle) and middle[n_tau].split(",", 1)[0] == w0:
            n_tau += 1
    blocks["tau"], blocks["induced_tau"] = middle[:n_tau], middle[n_tau:]
    for ch in ("Dmeson", "Bc", "Kmeson"):
        if own[ch] != blocks[ch]:
            raise ValueError(f"{flavor} {mass}: the {ch} block of the stored combined file is not its channel file")
    if blocks["induced_tau"] and not blocks["tau"]:
        raise ValueError(f"{flavor} {mass}: induced-tau rows without a prompt-tau block")
    return {ch: blocks[ch] for ch in ORDER}


def to_array(rows: list[str]) -> np.ndarray:
    if not rows:
        return np.empty((0, 5))
    return np.loadtxt(rows, delimiter=",", ndmin=2)
