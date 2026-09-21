#!/usr/bin/env python3
"""Generate full-branching BC10 decay templates with ROOT/Pythia8.

The selected exact SensCalc exclusive branching ratios determine the primary
decay. Pythia decays unstable daughters and
showers/hadronizes partonic modes, so the cached template contains stable final
particles and the GRENDEL reconstruction determines visibility. Three-body
primaries are generated in flat phase space and carry normalized weights from
the corresponding exact exported matrix elements. Pythia cannot hadronize an
isolated two-gluon colour singlet through this external-decay interface, so
``a -> gg`` is represented by an equal u/d/s light-quark jet mixture. The
``--gluon-surrogate`` option generates the pure-u, pure-d, and pure-s
alternatives used to propagate this implementation choice as a
decay-acceptance uncertainty.

Run under the same PyROOT + Pythia8 environment as the HNL FairShip templates:

    python -m grendel.models.alp.templates_pythia
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

from ...io.paths import ModelPaths
from ...io.vectors import format_mass_for_filename
from ..hnl.templates.fairship_decay import (
    PythiaCommandAdapter,
    _extract_underlying_pythia,
    _particle_charge,
    _require_root,
    _tpythia8_instance,
)
from .decay_matrix_elements import normalized_template_weights
from .decay_models import DECAY_MODEL_SPECS, DEFAULT_DECAY_MODEL, validate_decay_model
from .exclusive_decays import (
    DECAY_PRODUCTS_PDG,
    UNCLASSIFIED_CHANNEL_ID,
    ctau_at_reference_coupling,
    exclusive_branching_weights,
)
from .mass_grid import ALP_MASS_GRID

ALP_PDG = 9900015
RESONANCE_WINDOWS = ((0.125, 0.140), (0.538, 0.555), (0.940, 0.974))
GLUON_CHANNEL_ID = "channel_018"

_PRODUCTS_TO_CHANNEL = {
    tuple(products): channel_id
    for channel_id, products in DECAY_PRODUCTS_PDG.items()
    if channel_id not in {"channel_024", "channel_029"}
}


def _excluded_resonance(mass_gev: float) -> bool:
    return any(lower < mass_gev < upper for lower, upper in RESONANCE_WINDOWS)


def _mass_label(mass_gev: float) -> str:
    return f"{mass_gev:.3f}".replace(".", "p")


def _masses_from_csv(path: Path) -> list[float]:
    with path.expanduser().resolve().open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or "mass_GeV" not in reader.fieldnames:
            raise ValueError(f"mass-grid file has no mass_GeV column: {path}")
        masses = [float(row["mass_GeV"]) for row in reader]
    values = np.asarray(masses, dtype=float)
    if len(values) == 0 or np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("mass grid must contain finite positive values")
    if np.any(np.diff(values) <= 0.0):
        raise ValueError("mass grid must be strictly increasing with no duplicates")
    return masses


def _validate_resumed_template(
    path: Path,
    n_templates: int,
    gluon_surrogate: str,
    decay_model: str,
) -> None:
    """Refuse to reuse a template generated for a different model/config."""
    with np.load(path) as bundle:
        observed_model = (
            str(bundle["decay_model"])
            if "decay_model" in bundle.files
            else DEFAULT_DECAY_MODEL
        )
        observed_count = int(bundle["n_templates"])
        observed_surrogate = str(bundle["gluon_surrogate"])
    expected = (decay_model, int(n_templates), gluon_surrogate)
    observed = (observed_model, observed_count, observed_surrogate)
    if observed != expected:
        raise RuntimeError(
            f"resume template {path} has model/count/surrogate {observed}; "
            f"expected {expected}"
        )


class AlpPythiaBackend:
    """One initialized Pythia instance for an ALP mass and BR mixture."""

    def __init__(
        self,
        mass_gev: float,
        seed: int,
        gluon_surrogate: str = "uds",
        decay_model: str = DEFAULT_DECAY_MODEL,
    ):
        self.mass_gev = float(mass_gev)
        if gluon_surrogate not in {"uds", "u", "d", "s"}:
            raise ValueError(f"unknown gluon surrogate {gluon_surrogate!r}")
        self.gluon_surrogate = gluon_surrogate
        self.decay_model = validate_decay_model(decay_model)
        self.root = _require_root()
        self.tp8 = _tpythia8_instance(self.root)
        if self.tp8 is None:
            raise RuntimeError("ROOT.TPythia8 is required")
        self.adapter = PythiaCommandAdapter(self.tp8)
        self.engine = _extract_underlying_pythia(self.tp8)
        self.weights = exclusive_branching_weights(
            self.mass_gev, self.decay_model
        )
        self._configure(seed)

    def _read(self, command: str) -> None:
        if hasattr(self.engine, "readString"):
            self.engine.readString(command)
        else:
            self.tp8.ReadString(command)

    def _configure(self, seed: int) -> None:
        for command in (
            "ProcessLevel:all = off",
            "Print:quiet = on",
            "Init:showChangedSettings = off",
            "Init:showChangedParticleData = off",
            "Next:numberShowInfo = 0",
            "Next:numberShowProcess = 0",
            "Next:numberShowEvent = 0",
            "Random:setSeed = on",
            f"Random:seed = {seed}",
        ):
            self._read(command)
        self.adapter.SetParameters(
            f"{ALP_PDG}:new = a a 1 0 0 {self.mass_gev:.12g} 0 0 0 0 0 1 0 1 0"
        )
        self.adapter.SetParameters(f"{ALP_PDG}:isResonance = false")
        self.adapter.SetParameters(f"{ALP_PDG}:onMode = off")
        for channel_id, branching in self.weights.items():
            if channel_id == GLUON_CHANNEL_ID:
                # The external-decay interface does not create a fragmentable
                # colour topology for a bare gg pair. A light-flavour mixture
                # retains a conservative low-mass jet multiplicity model.
                flavors = {
                    "uds": (1, 2, 3),
                    "u": (2,),
                    "d": (1,),
                    "s": (3,),
                }[self.gluon_surrogate]
                modes = tuple(
                    (branching / len(flavors), (q, -q)) for q in flavors
                )
            else:
                products = (
                    (22, 22)
                    if channel_id == UNCLASSIFIED_CHANNEL_ID
                    else DECAY_PRODUCTS_PDG[channel_id]
                )
                modes = ((branching, products),)
            for mode_branching, products in modes:
                codes = " ".join(str(code) for code in products)
                self.adapter.SetParameters(
                    f"{ALP_PDG}:addChannel = 1 {mode_branching:.16g} 0 {codes}"
                )
        self.adapter.SetParameters(f"{ALP_PDG}:mayDecay = on")
        if hasattr(self.engine, "init") and not self.engine.init():
            raise RuntimeError(f"Pythia initialization failed at {self.mass_gev:g} GeV")

    def sample(self) -> dict[str, np.ndarray]:
        event = self.engine.event
        event.reset()
        event.append(
            ALP_PDG, 1, 0, 0, 0.0, 0.0, 0.0,
            self.mass_gev, self.mass_gev, 0.0, 9.0,
        )
        if not self.engine.next():
            raise RuntimeError(f"Pythia decay failed at {self.mass_gev:g} GeV")
        parents = [
            index for index in range(event.size())
            if int(event[index].id()) == ALP_PDG
        ]
        if len(parents) != 1:
            raise RuntimeError("Pythia event does not contain exactly one ALP")
        parent = parents[0]
        daughter_first = int(event[parent].daughter1())
        daughter_last = int(event[parent].daughter2())
        primary = tuple(
            int(event[index].id())
            for index in range(daughter_first, daughter_last + 1)
        )
        channel_id = _PRODUCTS_TO_CHANNEL.get(primary, "unweighted")
        primary_e1 = float(event[daughter_first].e()) if len(primary) == 3 else np.nan
        primary_e3 = float(event[daughter_last].e()) if len(primary) == 3 else np.nan
        rows = []
        for index in range(event.size()):
            particle = event[index]
            ancestor = int(particle.mother1())
            while ancestor > 0 and ancestor != parent:
                ancestor = int(event[ancestor].mother1())
            is_leaf = int(particle.daughter1()) <= 0 and int(particle.daughter2()) <= 0
            if ancestor != parent or not is_leaf:
                continue
            pdg = int(particle.id())
            rows.append((
                pdg,
                float(particle.px()),
                float(particle.py()),
                float(particle.pz()),
                float(particle.e()),
                float(particle.m()),
                _particle_charge(self.root, pdg),
            ))
        if not rows:
            raise RuntimeError("Pythia returned no stable ALP descendants")
        data = np.asarray(rows, dtype=float)
        return {
            "pdg": data[:, 0].astype(np.int32),
            "px": data[:, 1],
            "py": data[:, 2],
            "pz": data[:, 3],
            "energy": data[:, 4],
            "mass": data[:, 5],
            "charge": data[:, 6],
            "stable": np.ones(len(data), dtype=bool),
            "primary_channel": channel_id,
            "primary_energy_1": primary_e1,
            "primary_energy_3": primary_e3,
        }


def _flatten(samples: list[dict[str, np.ndarray]]) -> dict[str, np.ndarray]:
    counts = np.asarray([len(sample["pdg"]) for sample in samples], dtype=np.int32)
    return {
        "daughter_counts": counts,
        **{
            key: np.concatenate([sample[key] for sample in samples])
            for key in ("pdg", "px", "py", "pz", "energy", "mass", "charge", "stable")
        },
    }


def generate_one(
    mass_gev: float,
    n_templates: int,
    out_dir: Path,
    seed: int,
    gluon_surrogate: str = "uds",
    decay_model: str = DEFAULT_DECAY_MODEL,
    resume: bool = False,
) -> Path | None:
    if _excluded_resonance(mass_gev):
        print(f"  m_a={mass_gev:.3f}: skip unsupported light-meson resonance")
        return None
    destination = out_dir / f"templates_{_mass_label(mass_gev)}.npz"
    if resume and destination.exists():
        _validate_resumed_template(
            destination, n_templates, gluon_surrogate, decay_model
        )
        print(f"  m_a={mass_gev:.3f}: already checkpointed -> {destination.name}")
        return destination
    backend = AlpPythiaBackend(
        mass_gev, seed, gluon_surrogate, decay_model
    )
    samples = [backend.sample() for _ in range(n_templates)]
    bundle = _flatten(samples)
    primary_channels = np.asarray(
        [sample["primary_channel"] for sample in samples], dtype="U32"
    )
    primary_energy_1 = np.asarray(
        [sample["primary_energy_1"] for sample in samples], dtype=float
    )
    primary_energy_3 = np.asarray(
        [sample["primary_energy_3"] for sample in samples], dtype=float
    )
    matrix_element_weight = normalized_template_weights(
        primary_channels,
        primary_energy_1,
        primary_energy_3,
        mass_gev,
        decay_model,
    )
    charged = []
    offset = 0
    for count in bundle["daughter_counts"]:
        end = offset + int(count)
        charged.append(int((np.abs(bundle["charge"][offset:end]) > 0.5).sum()))
        offset = end
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(
            stream,
            **bundle,
            mass_GeV=np.float64(mass_gev),
            ctau_m_u2eq1=np.float64(
                ctau_at_reference_coupling(mass_gev, decay_model)
            ),
            n_templates=np.int32(n_templates),
            seed=np.int64(seed),
            flavor=np.array("BC10"),
            includes_full_branching=np.bool_(True),
            primary_channel=primary_channels,
            primary_energy_1_GeV=primary_energy_1,
            primary_energy_3_GeV=primary_energy_3,
            matrix_element_weight=matrix_element_weight,
            decay_backend=np.array(
                f"Pythia8-{decay_model}-weighted-three-body-"
                f"gg-to-{gluon_surrogate}"
            ),
            decay_model=np.array(decay_model),
            gluon_surrogate=np.array(gluon_surrogate),
        )
    temporary.replace(destination)
    print(
        f"  m_a={mass_gev:.3f}: stable-track visible="
        f"{np.mean(np.asarray(charged) >= 2):.2%} -> {destination.name}",
        flush=True,
    )
    return destination


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mass_source = parser.add_mutually_exclusive_group()
    mass_source.add_argument("--mass", type=float, nargs="+", default=None)
    mass_source.add_argument(
        "--mass-grid-file",
        type=Path,
        default=None,
        help="CSV containing the exact strictly increasing mass_GeV list",
    )
    parser.add_argument("--n-templates", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--out", type=Path, default=None,
                        help="template directory (default: $GRENDEL_BC10_TEMPLATES_DIR)")
    parser.add_argument(
        "--decay-model",
        choices=tuple(DECAY_MODEL_SPECS),
        default=DEFAULT_DECAY_MODEL,
        help="pinned decay model used for widths, BRs, and matrix elements",
    )
    parser.add_argument(
        "--gluon-surrogate",
        choices=("uds", "u", "d", "s"),
        default="uds",
        help="partonic proxy used for a -> gg (default: equal u/d/s mixture)",
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="keep atomically completed per-mass template bundles",
    )
    args = parser.parse_args(argv)
    try:
        if args.mass is not None:
            masses = args.mass
        elif args.mass_grid_file is not None:
            masses = _masses_from_csv(args.mass_grid_file)
        else:
            masses = ALP_MASS_GRID
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    for mass in masses:
        generate_one(
            mass,
            args.n_templates,
            args.out,
            args.seed + int(round(mass * 1000)),
            args.gluon_surrogate,
            args.decay_model,
            args.resume,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
