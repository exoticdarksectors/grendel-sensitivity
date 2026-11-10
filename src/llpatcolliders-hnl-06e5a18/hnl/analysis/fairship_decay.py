"""
analysis/fairship_decay.py

FairShip-backed HNL lifetime and decay sampling for the GARGOYLE analysis.

Supports all three flavors (Ue, Umu, Utau) via the FLAVOR_TO_COUPLINGS mapping.
Uses vendored FairShip modules (hnl.py, readDecayTable.py, shipunit.py,
pythia8_conf_utils.py) instead of requiring a full FairShip installation.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
VENDORED_FAIRSHIP = PROJECT_ROOT / "vendored" / "fairship"

HNL_PDG = 9900015
DEFAULT_DECAY_CONFIG = Path(__file__).with_name("fairship_decay_selection.conf")
SPEED_OF_LIGHT = 299_792_458.0  # m / s

FLAVOR_TO_COUPLINGS = {
    "Ue":   [1.0, 0.0, 0.0],
    "Umu":  [0.0, 1.0, 0.0],
    "Utau": [0.0, 0.0, 1.0],
}

_ROOT = None
_ROOT_ERROR = None
_MODULE_CACHE = {}
_BACKEND_CACHE = {}


@dataclass(frozen=True)
class DecayTemplate:
    """Rest-frame HNL decay template."""

    pdg: np.ndarray
    px: np.ndarray
    py: np.ndarray
    pz: np.ndarray
    energy: np.ndarray
    mass: np.ndarray
    charge: np.ndarray
    stable: np.ndarray
    status: np.ndarray


def _ensure_fairship_paths():
    """Insert vendored fairship directory into sys.path and set FAIRSHIP env."""
    vendored_str = str(VENDORED_FAIRSHIP)
    os.environ.setdefault("FAIRSHIP", vendored_str)
    if vendored_str not in sys.path:
        sys.path.insert(0, vendored_str)


def _try_import_root():
    global _ROOT, _ROOT_ERROR

    if _ROOT is not None:
        return _ROOT

    try:
        import ROOT  # type: ignore

        _ROOT = ROOT
        return ROOT
    except Exception as first_error:  # pragma: no cover - depends on local ROOT
        root_config = shutil.which("root-config")
        if root_config:
            try:
                libdir = subprocess.run(
                    [root_config, "--libdir"],
                    check=True,
                    capture_output=True,
                    text=True,
                ).stdout.strip()
                if libdir and libdir not in sys.path:
                    sys.path.insert(0, libdir)
                import ROOT  # type: ignore

                _ROOT = ROOT
                return ROOT
            except Exception as second_error:  # pragma: no cover - env specific
                _ROOT_ERROR = second_error
        if _ROOT_ERROR is None:
            _ROOT_ERROR = first_error
        return None


def _require_root():
    global _ROOT

    root = _try_import_root()
    if root is None:
        py_version = "unknown"
        root_config = shutil.which("root-config")
        if root_config:
            try:
                py_version = subprocess.run(
                    [root_config, "--python-version"],
                    check=True,
                    capture_output=True,
                    text=True,
                ).stdout.strip()
            except Exception:
                py_version = "unknown"
        raise RuntimeError(
            "FairShip decay backend requires a working PyROOT + Pythia8 setup. "
            f"Could not import ROOT ({_ROOT_ERROR}). "
            f"This ROOT build expects Python {py_version}. "
            "Activate the llpatcolliders conda environment before running."
        ) from _ROOT_ERROR

    _ROOT = root
    _load_pythia_libraries(root)
    return root


def _load_pythia_libraries(root):
    for lib_name in ("libEG", "libPythia8", "libpythia8", "libEGPythia8"):
        try:  # pragma: no branch - library availability is environment specific
            root.gSystem.Load(lib_name)
        except Exception:
            continue


def _load_fairship_module(name):
    if name in _MODULE_CACHE:
        return _MODULE_CACHE[name]

    _ensure_fairship_paths()
    _require_root()
    module = __import__(name)
    _MODULE_CACHE[name] = module
    return module


def _hnl_module():
    return _load_fairship_module("hnl")


def _read_decay_table_module():
    return _load_fairship_module("readDecayTable")


def _pythia8_conf_utils_module():
    return _load_fairship_module("pythia8_conf_utils")


def _shipunit_module():
    return _load_fairship_module("shipunit")


def _couplings_for_flavor(flavor):
    """Return [U2e, U2mu, U2tau] coupling list for the given flavor."""
    if flavor in FLAVOR_TO_COUPLINGS:
        return FLAVOR_TO_COUPLINGS[flavor]
    raise ValueError(f"Unknown flavor: {flavor}. Use one of {list(FLAVOR_TO_COUPLINGS)}")


def _extract_underlying_pythia(tp8):
    if hasattr(tp8, "Pythia8"):
        return tp8.Pythia8()
    return tp8


class PythiaCommandAdapter:
    """Expose ``SetParameters(cmd)`` for FairShip's decay-table helpers."""

    def __init__(self, tp8):
        self._tp8 = tp8

    def _engine(self):
        return _extract_underlying_pythia(self._tp8)

    def SetParameters(self, cmd):
        engine = self._engine()
        if hasattr(engine, "readString"):
            return engine.readString(cmd)
        if hasattr(self._tp8, "ReadString"):
            return self._tp8.ReadString(cmd)
        raise AttributeError("Configured TPythia8 object cannot consume readString commands")

    def SetHNLId(self, pid):
        if hasattr(self._tp8, "SetHNLId"):
            return self._tp8.SetHNLId(pid)
        return None

    def List(self, pid):
        if hasattr(self._tp8, "List"):
            return self._tp8.List(pid)
        return None


def _tpythia8_instance(root):
    if not hasattr(root, "TPythia8"):
        return None

    try:
        instance = root.TPythia8.Instance()
    except Exception:
        instance = None
    if instance:
        return instance

    try:
        return root.TPythia8(False)
    except Exception:
        pass

    try:
        return root.TPythia8()
    except Exception:
        return None


def _particle_charge(root, pdg_code):
    particle = root.TDatabasePDG.Instance().GetParticle(int(pdg_code))
    if particle is None:
        return 0.0
    return float(particle.Charge()) / 3.0


def _particle_mass(root, pdg_code, energy, px, py, pz):
    particle = root.TDatabasePDG.Instance().GetParticle(int(pdg_code))
    if particle:
        return float(particle.Mass())

    mass2 = float(energy) ** 2 - (float(px) ** 2 + float(py) ** 2 + float(pz) ** 2)
    return float(np.sqrt(max(mass2, 0.0)))


def _register_hnl_in_root(root, hnl_instance, mass_gev):
    pdg = root.TDatabasePDG.Instance()
    if pdg.GetParticle(HNL_PDG):
        return

    units = _shipunit_module()
    ctau = hnl_instance.computeNLifetime(system="FairShip") * units.c_light * units.cm
    gamma = units.hbarc / float(ctau)
    _pythia8_conf_utils_module().addHNLtoROOT(pid=HNL_PDG, m=float(mass_gev), g=float(gamma))


def _build_hnl_instance(mass_gev, flavor="Umu"):
    """Build an HNL instance with the given mass and flavor couplings."""
    hnl = _hnl_module()
    couplings = _couplings_for_flavor(flavor)
    return hnl.HNL(mass_gev, couplings)


def _charge_conjugate_codes(root, codes):
    conjugates = []
    pdg = root.TDatabasePDG.Instance()
    for code in codes:
        anti = pdg.GetParticle(-int(code))
        conjugates.append(int(-code) if anti else int(code))
    return conjugates


def _add_hnl_decay_channels(adapter, hnl_instance, conffile, root):
    read_decay_table = _read_decay_table_module()
    allowed = hnl_instance.allowedChannels()
    wanted = read_decay_table.load(conffile=str(conffile), verbose=False)

    for decay, allowed_flag in allowed.items():
        if decay not in wanted:
            raise RuntimeError(f"HNL decay channel missing from decay config: {decay}")
        if allowed_flag != "yes" or wanted[decay] != "yes":
            continue

        particles = [particle for particle in decay.replace("->", " ").split()]
        children = particles[1:]
        child_codes = [int(read_decay_table.PDGcode(child)) for child in children]
        # Majorana HNL: each decay channel splits equally between the particle
        # and charge-conjugate final state (N = N-bar, so both contribute).
        # HNLCalc gives total BR; divide by 2 for each CP-conjugate pair.
        # For Dirac HNL this would need to be modified.
        br = float(hnl_instance.findBranchingRatio(decay)) / 2.0
        codes = " ".join(str(code) for code in child_codes)
        adapter.SetParameters(f"{HNL_PDG}:addChannel = 1 {br:.12} 0 {codes}")

        conj_codes = " ".join(str(code) for code in _charge_conjugate_codes(root, child_codes))
        adapter.SetParameters(f"{HNL_PDG}:addChannel = 1 {br:.12} 0 {conj_codes}")


def ctau_u2eq1(mass_GeV, flavor="Umu"):
    """Return the proper decay length at ``U^2 = 1`` in metres."""
    hnl_instance = _build_hnl_instance(mass_GeV, flavor=flavor)
    return SPEED_OF_LIGHT * float(hnl_instance.computeNLifetime(system="SI"))


class FairShipDecayBackend:
    """Mass-point scoped FairShip decay backend."""

    def __init__(self, mass_GeV, flavor="Umu", decay_config=DEFAULT_DECAY_CONFIG, random_seed=None):
        self.mass_GeV = float(mass_GeV)
        self.flavor = flavor
        self.decay_config = Path(decay_config)
        self.random_seed = None if random_seed is None else int(random_seed)

        self.root = _require_root()
        self.hnl = _build_hnl_instance(self.mass_GeV, flavor=self.flavor)
        self.ctau_m = ctau_u2eq1(self.mass_GeV, flavor=self.flavor)
        self._decayer = None
        self._tp8 = _tpythia8_instance(self.root)
        if self._tp8 is None:
            raise RuntimeError(
                "PyROOT is available, but ROOT.TPythia8 is not. "
                "This analysis requires a ROOT build with Pythia8 support."
            )

        self._adapter = PythiaCommandAdapter(self._tp8)
        self._engine = _extract_underlying_pythia(self._tp8)

        self._configure_pythia()

    def _configure_pythia(self):
        quiet_commands = (
            "ProcessLevel:all = off",
            "Print:quiet = on",
            "Init:showChangedSettings = off",
            "Init:showChangedParticleData = off",
            "Next:numberShowInfo = 0",
            "Next:numberShowProcess = 0",
            "Next:numberShowEvent = 0",
        )
        for command in quiet_commands:
            if hasattr(self._tp8, "ReadString"):
                self._tp8.ReadString(command)
            elif hasattr(self._engine, "readString"):
                self._engine.readString(command)

        if self.random_seed is not None and hasattr(self._engine, "readString"):
            self._engine.readString("Random:setSeed = on")
            self._engine.readString(f"Random:seed = {self.random_seed}")

        ctau_mm = self.ctau_m * 1.0e3
        self._adapter.SetParameters(
            f"{HNL_PDG}:new = N2 N2 2 0 0 {self.mass_GeV:.12} 0.0 0.0 0.0 {ctau_mm:.12} 0 1 0 1 0"
        )
        self._adapter.SetParameters(f"{HNL_PDG}:isResonance = false")
        self._adapter.SetParameters(f"{HNL_PDG}:onMode = off")
        _add_hnl_decay_channels(self._adapter, self.hnl, self.decay_config, self.root)
        self._adapter.SetParameters(f"{HNL_PDG}:mayDecay = on")
        self._adapter.SetHNLId(HNL_PDG)
        _register_hnl_in_root(self.root, self.hnl, self.mass_GeV)
        if hasattr(self._engine, "init"):
            self._engine.init()

    def _sample_with_decayer(self):
        if self._decayer is None:
            if not hasattr(self.root, "TPythia8Decayer"):
                raise RuntimeError("ROOT.TPythia8Decayer is not available in this ROOT build")
            self._decayer = self.root.TPythia8Decayer()
            self._decayer.Init()

        particles = self.root.TClonesArray("TParticle")
        parent = self.root.TLorentzVector(0.0, 0.0, 0.0, self.mass_GeV)
        self._decayer.Decay(HNL_PDG, parent)
        n_particles = int(self._decayer.ImportParticles(particles))
        return self._particles_to_template(particles, n_particles)

    def _sample_with_raw_pythia(self):  # pragma: no cover - depends on TPythia8 bindings
        event = self._engine.event
        event.reset()
        event.append(HNL_PDG, 1, 0, 0, 0.0, 0.0, 0.0, self.mass_GeV, self.mass_GeV, 0.0, 9.0)
        self._engine.next()

        hnl_indices = [idx for idx in range(event.size()) if abs(int(event[idx].id())) == HNL_PDG]
        if not hnl_indices:
            raise RuntimeError(f"Pythia did not produce an HNL decay at m_N={self.mass_GeV:.3f} GeV")
        hnl_idx = hnl_indices[0]

        pdg, px, py, pz, energy, mass, charge, stable, status = ([] for _ in range(9))
        for idx in range(event.size()):
            particle = event[idx]
            ancestor = int(particle.mother1())
            while ancestor > 0 and ancestor != hnl_idx:
                ancestor = int(event[ancestor].mother1())
            daughter1 = int(particle.daughter1())
            daughter2 = int(particle.daughter2())
            is_leaf = daughter1 <= 0 and daughter2 <= 0
            if ancestor != hnl_idx or not is_leaf:
                continue

            pdg_code = int(particle.id())
            pdg.append(pdg_code)
            px.append(float(particle.px()))
            py.append(float(particle.py()))
            pz.append(float(particle.pz()))
            energy.append(float(particle.e()))
            mass.append(float(particle.m()))
            charge.append(_particle_charge(self.root, pdg_code))
            stable.append(True)
            status.append(1)

        return DecayTemplate(
            pdg=np.asarray(pdg, dtype=np.int32),
            px=np.asarray(px, dtype=np.float64),
            py=np.asarray(py, dtype=np.float64),
            pz=np.asarray(pz, dtype=np.float64),
            energy=np.asarray(energy, dtype=np.float64),
            mass=np.asarray(mass, dtype=np.float64),
            charge=np.asarray(charge, dtype=np.float64),
            stable=np.asarray(stable, dtype=bool),
            status=np.asarray(status, dtype=np.int32),
        )

    def _particles_to_template(self, particles, n_particles):
        pdg, px, py, pz, energy, mass, charge, stable, status = ([] for _ in range(9))
        for idx in range(n_particles):
            particle = particles.At(idx)
            if particle is None:
                continue

            pdg_code = int(particle.GetPdgCode())
            if abs(pdg_code) == HNL_PDG:
                continue

            status_code = int(particle.GetStatusCode())
            px_i = float(particle.Px())
            py_i = float(particle.Py())
            pz_i = float(particle.Pz())
            energy_i = float(particle.Energy())
            pdg.append(pdg_code)
            px.append(px_i)
            py.append(py_i)
            pz.append(pz_i)
            energy.append(energy_i)
            mass.append(_particle_mass(self.root, pdg_code, energy_i, px_i, py_i, pz_i))
            charge.append(_particle_charge(self.root, pdg_code))
            stable.append(status_code == 1)
            status.append(status_code)

        return DecayTemplate(
            pdg=np.asarray(pdg, dtype=np.int32),
            px=np.asarray(px, dtype=np.float64),
            py=np.asarray(py, dtype=np.float64),
            pz=np.asarray(pz, dtype=np.float64),
            energy=np.asarray(energy, dtype=np.float64),
            mass=np.asarray(mass, dtype=np.float64),
            charge=np.asarray(charge, dtype=np.float64),
            stable=np.asarray(stable, dtype=bool),
            status=np.asarray(status, dtype=np.int32),
        )

    def sample_rest_frame_decay(self):
        template = self._sample_with_raw_pythia()
        if len(template.pdg) == 0:
            raise RuntimeError(f"FairShip decay sampling returned no daughters at m_N={self.mass_GeV:.3f} GeV")
        return template

    def sample_rest_frame_decays(self, n_templates):
        return [self.sample_rest_frame_decay() for _ in range(int(n_templates))]


def _backend_cache_key(mass_GeV, flavor, decay_config, random_seed):
    return (round(float(mass_GeV), 12), flavor, str(Path(decay_config).resolve()), random_seed)


_MAX_BACKEND_CACHE = 8


def get_decay_backend(mass_GeV, flavor="Umu", decay_config=DEFAULT_DECAY_CONFIG, random_seed=None):
    if random_seed is not None:
        return FairShipDecayBackend(
            mass_GeV=mass_GeV,
            flavor=flavor,
            decay_config=decay_config,
            random_seed=random_seed,
        )

    key = _backend_cache_key(mass_GeV, flavor, decay_config, random_seed)
    if key not in _BACKEND_CACHE:
        if len(_BACKEND_CACHE) >= _MAX_BACKEND_CACHE:
            _BACKEND_CACHE.pop(next(iter(_BACKEND_CACHE)))
        _BACKEND_CACHE[key] = FairShipDecayBackend(
            mass_GeV=mass_GeV,
            flavor=flavor,
            decay_config=decay_config,
            random_seed=random_seed,
        )
    return _BACKEND_CACHE[key]


def sample_rest_frame_decays(mass_GeV, n_templates, flavor="Umu",
                             decay_config=DEFAULT_DECAY_CONFIG, random_seed=None):
    """Return a list of FairShip-sampled rest-frame HNL decays."""
    return get_decay_backend(
        mass_GeV, flavor=flavor, decay_config=decay_config, random_seed=random_seed,
    ).sample_rest_frame_decays(n_templates)


def boost_decay_to_lab(parent_p4, template):
    """
    Boost a rest-frame decay template into the supplied lab-frame parent 4-vector.

    Uses a general Lorentz boost along the actual parent 3-momentum direction
    (NOT a z-only boost).

    Parameters
    ----------
    parent_p4 : array-like of shape (4,)
        ``(E, px, py, pz)`` in GeV.
    template : DecayTemplate

    Returns
    -------
    dict
        Daughter 4-vectors and metadata in the lab frame.
    """
    parent_p4 = np.asarray(parent_p4, dtype=np.float64)
    if parent_p4.shape != (4,):
        raise ValueError("parent_p4 must have shape (4,) with (E, px, py, pz)")

    energy = parent_p4[0]
    momentum = parent_p4[1:]
    if energy <= 0:
        raise ValueError("parent energy must be positive")

    beta = momentum / energy
    beta2 = float(np.dot(beta, beta))
    if beta2 >= 1.0:
        beta2 = np.nextafter(1.0, 0.0)

    rest_momenta = np.column_stack([template.px, template.py, template.pz])
    rest_energies = template.energy
    if beta2 <= 0.0:
        lab_momenta = rest_momenta.copy()
        lab_energies = rest_energies.copy()
    else:
        gamma = 1.0 / np.sqrt(1.0 - beta2)
        bp = rest_momenta @ beta
        gamma2 = (gamma - 1.0) / beta2
        lab_momenta = rest_momenta + (gamma2 * bp + gamma * rest_energies)[:, None] * beta[None, :]
        lab_energies = gamma * (rest_energies + bp)

    return {
        "pdg": template.pdg,
        "px": lab_momenta[:, 0],
        "py": lab_momenta[:, 1],
        "pz": lab_momenta[:, 2],
        "energy": lab_energies,
        "mass": template.mass,
        "charge": template.charge,
        "stable": template.stable,
        "status": template.status,
    }
