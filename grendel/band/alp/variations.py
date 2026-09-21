"""The variations of the BC10 uncertainty campaign.

* the coherent FONLL bottom ensemble: central, six scale grids, 100 NNPDF
  replicas, two bottom-mass grids -- each with a fresh 600k production;
* three ``a -> gg`` hadronisation surrogates (pure u, d, s light-quark
  jets instead of the equal mixture) on the central production vectors;
* the flavour-violating coefficient ``C_bs`` at +/-20% amplitude
  (production rates scaled by 0.8^2 and 1.2^2) with fresh production;
* two same-physics central repeats with independent production and
  reconstruction seeds, the numerical controls.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

from ...io.atomic import sha256_file
from ..fonll_variations import assert_counts, bottom_records

CBS_SCHEME_RELATIVE_AMPLITUDE = 0.20
EXPECTED_FONLL = {"central": 1, "scale": 6, "pdf": 100, "mass": 2}


def stable_seed(name: str) -> int:
    """Stable independent NumPy seed for one non-central production run."""
    digest = hashlib.sha256(f"GRENDEL-BC10-FONLL:{name}".encode()).digest()
    return int.from_bytes(digest[:4], "big") or 1


def discover_fonll_variations(grid_dir: Path) -> list[dict]:
    """The complete bottom-grid ensemble from its pinned manifest."""
    records = bottom_records(grid_dir, mass_naming="tag")
    assert_counts(records, EXPECTED_FONLL)
    variations = []
    for r in records:
        variations.append({
            "name": r["name"],
            "axis": "mb" if r["axis"] == "mass" else r["axis"],
            "campaign_axis": "fonll",
            "grid_path": str(r["path"].resolve()),
            "grid_sha256": r["sha256"],
            "muR": r["muR"], "muF": r["muF"],
            "pdf_member": r["lhapdf_member"],
            "mb_GeV": r["heavy_quark_mass_GeV"],
            "production_seed": 42 if r["axis"] == "central" else stable_seed(r["tag"]),
            "reco_seed_offset": 0,
            "cbs_amplitude_scale": 1.0,
            "template_variant": "central",
            "production_mode": "fresh_600k",
        })
    return variations


def auxiliary_variations(central_grid: Path) -> list[dict]:
    grid_path = str(Path(central_grid).resolve())
    grid_sha = sha256_file(Path(central_grid))
    output = []
    for flavor in ("u", "d", "s"):
        output.append({"name": f"gg_{flavor}", "axis": "decay_gg", "campaign_axis": "decay_gg",
                       "grid_path": grid_path, "grid_sha256": grid_sha, "production_seed": 42,
                       "reco_seed_offset": 0, "cbs_amplitude_scale": 1.0,
                       "template_variant": f"gg_{flavor}",
                       "production_mode": "central_vectors_exact_reuse"})
    for direction, scale in (("down", 1.0 - CBS_SCHEME_RELATIVE_AMPLITUDE),
                             ("up", 1.0 + CBS_SCHEME_RELATIVE_AMPLITUDE)):
        output.append({"name": f"cbs_{direction}", "axis": "cbs", "campaign_axis": "cbs",
                       "grid_path": grid_path, "grid_sha256": grid_sha, "production_seed": 42,
                       "reco_seed_offset": 0, "cbs_amplitude_scale": scale,
                       "template_variant": "central", "production_mode": "fresh_600k"})
    return output


def numerical_control_variations(central_grid: Path) -> list[dict]:
    """Same-physics central repeats used only to measure numerical spread."""
    grid_path = str(Path(central_grid).resolve())
    grid_sha = sha256_file(Path(central_grid))
    output = []
    for index in (1, 2):
        name = f"central_repeat_{index}"
        output.append({"name": name, "axis": "numerical_control", "campaign_axis": "numerical_control",
                       "grid_path": grid_path, "grid_sha256": grid_sha,
                       "production_seed": stable_seed(f"{name}:production"),
                       "reco_seed_offset": stable_seed(f"{name}:reconstruction"),
                       "cbs_amplitude_scale": 1.0, "template_variant": "central",
                       "production_mode": "fresh_600k"})
    return output


def all_variations(grid_dir: Path) -> list[dict]:
    fonll = discover_fonll_variations(grid_dir)
    central_grid = Path(fonll[0]["grid_path"])
    return [*fonll, *auxiliary_variations(central_grid), *numerical_control_variations(central_grid)]
