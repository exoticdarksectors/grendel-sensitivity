"""Generate the FONLL meson grids: central, or the variation campaign.

    python -m grendel.production.fonll_grids.generate --pdf nlo --quark bottom --grid-workers 4
    python -m grendel.production.fonll_grids.generate --campaign --max-parallel 6 --compress-logs

Each grid runs the patched FONLL ``fonllgridlha`` (the heavy-quark grid,
optionally split along rapidity into ``--grid-workers`` independent
processes and re-stitched) and then ``fragmfonll`` (the fragmentation
convolution on the 100 x 100 (pT, y) meson grid). Bottom is the public
FONLL B-hadron default (Kartvelishvili, alpha = 24.2). Charm is the public
D0 convention: the direct BCFY pseudoscalar stream plus a collinear
D* -> D0 feeddown of the BCFY vector stream, combined with a single weight
fitted once to nine public FONLL CTEQ6.6 D0 points and cached.

``--campaign`` enumerates the variation space (7-point scale, NNPDF4.0
replica members, heavy-quark mass) and writes ``variation_manifest.json``
with a SHA-256 per grid; the band campaigns read that manifest.

Locations: ``--fonll`` (``GRENDEL_FONLL_DIR``) is the FONLL tree with the
executables built under ``Linux/``; ``--out-dir`` (``GRENDEL_FONLL_GRID_DIR``)
receives the grids; run and log directories default beside it.
``GRENDEL_LHAPDF_DATA`` is exported as ``LHAPDF_DATA_PATH`` when set.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from ...io.paths import fonll_dir, fonll_grid_dir

PDFS = {
    "cteq66": ("CTEQ6.6", 10550, "cteq66"),
    "nlo": ("NNPDF40_nlo_as_01180", 331700),
    "nlo_as_01170": ("NNPDF40_nlo_as_01170", 333900, "nnpdf40_nlo_as_01170"),
    "nlo_as_01190": ("NNPDF40_nlo_as_01190", 334100, "nnpdf40_nlo_as_01190"),
    "nnlo": ("NNPDF40_nnlo_as_01180", 331100),
}

QUARKS = {
    "bottom": {
        "mass": 4.75,
        "frag_mode": 2,
        "frag_param": 24.2,
        "frag_param_name": "alpha",
        "fragmentation": "Kartvelishvili, normalized (1-z)*z^alpha; public FONLL N=5 central default for B hadrons",
        "final_state_detail": "B hadron",
    },
    "charm": {
        "mass": 1.50,
        "frag_mode": 5,
        "frag_param": 0.1,
        "frag_param_name": "r",
        "fragmentation": "BCFY pseudoscalar plus calibrated D* feeddown; public FONLL D0 convention",
        "final_state_detail": "D0 with public FONLL D* feeddown convention",
        "feeddown": "public_d0",
    },
}

# The meson grid resolution the samplers were built for; changing the PDF
# set must not silently alter it.
Y_VALUES = [round(-3.0 + 6.0 * i / 99.0, 10) for i in range(100)]
PT_VALUES = [round(50.0 * i / 99.0, 10) for i in range(100)]

# The internal quark grid is wider than the meson grid because
# fragmentation samples quark pT above the observed meson pT.
FONLL_GRID_PTMAX = 200.0
FONLL_GRID_NPT = 80
FRAG_FRAME = 1
FRAG_MESON_MASS = -1.0
MATCHING_C = 5.0

# Standard 7-point scale variation as (muR/mu0, muF/mu0), mu0 = sqrt(m^2 + pT^2),
# without the two antipodal extremes. FONLL's grid driver reads the factors
# in the order (ffact = muF, fren = muR): misc1/fonllgrid.f `read(*,*) ffact,fren`,
# main/fonll0.f (ffact -> xmu_fact, fren -> xmu_ren).
SCALE_POINTS_7 = [(1.0, 1.0), (2.0, 2.0), (0.5, 0.5), (2.0, 1.0), (1.0, 2.0), (0.5, 1.0), (1.0, 0.5)]

# Heavy-quark-mass variations (GeV), FONLL benchmark conventions:
# m_b = 4.75 +/- 0.25, m_c = 1.5 -/+ 0.2.
MASS_VARIATIONS = {"bottom": [4.50, 5.00], "charm": [1.30, 1.70]}

# Public FONLL v1.3.2 CTEQ6.6 D0 references at 14 TeV, y = 0, queried from
# the public form on 2026-05-29. They calibrate only the charm D* feeddown
# weight; the generated grid still uses the requested LHAPDF set.
PUBLIC_CHARM_D0_CTEQ66_Y0 = {1.0: 2.5497e8, 5.0: 3.5221e7, 8.0: 6.8479e6, 15.0: 5.0952e5, 22.0: 9.0635e4,
                             29.0: 2.4923e4, 36.0: 8.8653e3, 43.0: 3.7354e3, 50.0: 1.7742e3}

# PDG masses (GeV) and branching fractions of the collinear two-body
# D* feeddown; the overall normalisation is fixed by the calibration above.
CHARM_D0_FEEDDOWN_CHANNELS = [
    {"name": "D*0 -> D0 pi0", "branching_fraction": 0.647, "parent_mass": 2.00685, "daughter_mass": 1.86484,
     "spectator_mass": 0.134977},
    {"name": "D*0 -> D0 gamma", "branching_fraction": 0.353, "parent_mass": 2.00685, "daughter_mass": 1.86484,
     "spectator_mass": 0.0},
    {"name": "D*+ -> D0 pi+", "branching_fraction": 0.677, "parent_mass": 2.01026, "daughter_mass": 1.86484,
     "spectator_mass": 0.139570},
]

_COMPRESS_LOGS = False


@dataclass(frozen=True)
class Workspace:
    """Where the FONLL executables are and where a run writes."""
    fonll: Path
    out: Path
    run: Path
    logs: Path

    @classmethod
    def resolve(cls, fonll=None, out=None, run=None, logs=None) -> "Workspace":
        fonll = Path(fonll) if fonll else fonll_dir()
        out = Path(out) if out else fonll_grid_dir()
        return cls(fonll=fonll.expanduser(), out=out.expanduser(),
                   run=Path(run).expanduser() if run else out.parent / "run",
                   logs=Path(logs).expanduser() if logs else out.parent / "logs")

    @property
    def linux(self) -> Path:
        return self.fonll / "Linux"

    def executable(self, name: str) -> Path:
        path = self.linux / name
        if not path.is_file():
            raise FileNotFoundError(f"{path} not found; build FONLL first "
                                    f"(python -m grendel.production.fonll_grids.install --build)")
        return path

    def env(self) -> dict[str, str]:
        env = os.environ.copy()
        data = os.environ.get("GRENDEL_LHAPDF_DATA")
        if data:
            env["LHAPDF_DATA_PATH"] = data
        return env


@dataclass(frozen=True)
class Variation:
    """One point in the (scale, PDF member, heavy-quark mass) space.

    ``muR``/``muF`` are relative to mu0; ``pdf_member`` is added to the set's
    base LHAPDF id; ``mass`` is the heavy-quark mass actually used.
    """
    kind: str   # "central" | "scale" | "pdf" | "mass"
    tag: str    # filesystem-safe label used in filenames and the manifest
    muR: float
    muF: float
    pdf_member: int
    mass: float


def _fmt_float_tag(value: float) -> str:
    """2.0 -> '2', 0.5 -> '0p5', 4.75 -> '4p75'."""
    return f"{value:g}".replace(".", "p").replace("-", "m")


def central_variation(quark: str) -> Variation:
    return Variation("central", "central", 1.0, 1.0, 0, float(QUARKS[quark]["mass"]))


def enumerate_variations(quark: str, pdf_members: Iterable[int], include_scale: bool,
                         include_mass: bool) -> list[Variation]:
    """The central point once (it serves as scale (1,1), member 0 and the
    central mass), then each axis's non-central points."""
    default_mass = float(QUARKS[quark]["mass"])
    variations = [central_variation(quark)]
    if include_scale:
        for muR, muF in SCALE_POINTS_7:
            if (muR, muF) == (1.0, 1.0):
                continue
            variations.append(Variation("scale", f"scale_muR{_fmt_float_tag(muR)}_muF{_fmt_float_tag(muF)}",
                                        muR, muF, 0, default_mass))
    for member in pdf_members:
        if member == 0:
            continue
        variations.append(Variation("pdf", f"pdf_{member:04d}", 1.0, 1.0, int(member), default_mass))
    if include_mass:
        for m in MASS_VARIATIONS[quark]:
            variations.append(Variation("mass", f"mass_{_fmt_float_tag(m)}", 1.0, 1.0, 0, float(m)))
    return variations


def point_key(pt: float, y: float) -> tuple[float, float]:
    # fragmfonll perturbs dense-grid coordinates at its printed precision;
    # five decimals stay far below the node spacing and keep every distinct
    # charm feeddown request apart.
    return (round(pt, 5), round(y, 5))


def daughter_energy_fraction(parent_mass: float, daughter_mass: float, spectator_mass: float) -> float:
    return (parent_mass ** 2 + daughter_mass ** 2 - spectator_mass ** 2) / (2.0 * parent_mass ** 2)


def charm_d0_feeddown_channels() -> list[dict[str, float | str]]:
    channels = []
    for channel in CHARM_D0_FEEDDOWN_CHANNELS:
        enriched = dict(channel)
        enriched["z_collinear"] = daughter_energy_fraction(
            float(channel["parent_mass"]), float(channel["daughter_mass"]), float(channel["spectator_mass"]))
        channels.append(enriched)
    return channels


def run_command(cmd: list[str], cwd: Path, stdin: str, log_path: Path, env: dict[str, str]) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        proc = subprocess.run(cmd, cwd=cwd, input=stdin, text=True, stdout=log, stderr=subprocess.STDOUT,
                              env=env, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"command failed with exit {proc.returncode}: {' '.join(cmd)}; see {log_path}")
    if _COMPRESS_LOGS:
        # Failing logs stay uncompressed (the error above points at the plain .log).
        gz_path = log_path.with_suffix(log_path.suffix + ".gz")
        with log_path.open("rb") as src, gzip.open(gz_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
        log_path.unlink()


def grid_input(prefix: str, pdf_id: int, mass: float, ffact: float = 1.0, fren: float = 1.0,
               y_values: Iterable[float] = Y_VALUES, ptmax: float = FONLL_GRID_PTMAX, npt: int = FONLL_GRID_NPT) -> str:
    """The ``fonllgridlha`` input card."""
    lines = [prefix, f" 1 7000. 0 0 {pdf_id}", f" 1 7000. 0 0 {pdf_id}", f" {mass:.8g}", " -1.",
             f" {ffact:.8g} {fren:.8g}"]      # (ffact = muF, fren = muR)
    lines.extend(f" {y:.8g}" for y in y_values)
    lines.extend([" 10000", " 0", f" {ptmax:.8g} {npt}", " 1"])
    return "\n".join(lines) + "\n"


def _contiguous_chunks(values: list[float], n_chunks: int) -> list[list[float]]:
    n_chunks = min(n_chunks, len(values))
    size, remainder = divmod(len(values), n_chunks)
    chunks, start = [], 0
    for index in range(n_chunks):
        stop = start + size + (1 if index < remainder else 0)
        chunks.append(values[start:stop])
        start = stop
    return chunks


def build_quark_grid(ws: Workspace, run_dir: Path, grid_file: Path, prefix: str, pdf_id: int, mass: float,
                     pdf_key: str, quark: str, grid_workers: int, ffact: float = 1.0, fren: float = 1.0,
                     log_label: str = "") -> None:
    env = ws.env()
    exe = str(ws.executable("fonllgridlha"))
    suffix = f"_{log_label}" if log_label else ""
    if grid_workers == 1:
        run_command([exe], run_dir, grid_input(prefix, pdf_id, mass, ffact=ffact, fren=fren),
                    ws.logs / f"fonllgrid_{pdf_key}_{quark}{suffix}.log", env)
        return
    chunk_root = run_dir / "grid_chunks"
    if chunk_root.exists():
        shutil.rmtree(chunk_root)
    chunk_root.mkdir()
    chunks = _contiguous_chunks(Y_VALUES, grid_workers)

    def run_chunk(index: int, y_values: list[float]) -> Path:
        chunk_dir = chunk_root / f"{index:02d}"
        chunk_dir.mkdir()
        run_command([exe], chunk_dir, grid_input(prefix, pdf_id, mass, ffact=ffact, fren=fren, y_values=y_values),
                    ws.logs / f"fonllgrid_{pdf_key}_{quark}{suffix}_chunk{index:02d}.log", env)
        chunk_file = chunk_dir / f"{prefix}.out"
        if not chunk_file.exists():
            raise RuntimeError(f"expected chunk grid file not found: {chunk_file}")
        return chunk_file

    with ThreadPoolExecutor(max_workers=len(chunks)) as executor:
        chunk_files = list(executor.map(run_chunk, range(len(chunks)), chunks))

    combined_lines: list[str] = []
    header = None
    for chunk_file in chunk_files:
        lines = chunk_file.read_text().splitlines()
        if not lines:
            raise RuntimeError(f"empty chunk grid file: {chunk_file}")
        if header is None:
            header = lines[0]
            combined_lines.append(header)
        elif lines[0] != header:
            raise RuntimeError(f"inconsistent raw-grid header in {chunk_file}")
        combined_lines.extend(lines[1:])
    expected_lines = 1 + len(Y_VALUES) * FONLL_GRID_NPT
    if len(combined_lines) != expected_lines:
        raise RuntimeError(f"merged raw grid has {len(combined_lines)} lines, expected {expected_lines}")
    grid_file.write_text("\n".join(combined_lines) + "\n")


def observed_points() -> list[tuple[float, float]]:
    return [(pt, y) for pt in PT_VALUES if pt > 0.0 for y in Y_VALUES]


def charm_feeddown_parent_points(points: Iterable[tuple[float, float]]) -> list[tuple[float, float]]:
    parent_points = set()
    for pt, y in points:
        for channel in charm_d0_feeddown_channels():
            parent_points.add(point_key(pt / float(channel["z_collinear"]), y))
    return sorted(parent_points)


def frag_input_points(grid_file_name: str, frag_mode: int, frag_param: float,
                      points: Iterable[tuple[float, float]]) -> str:
    """The ``fragmfonll`` input card for a list of (pT, y) evaluation points."""
    lines = ["0", "0", f"{MATCHING_C:.8g}", grid_file_name, "", "0", str(frag_mode), str(FRAG_FRAME),
             f"{FRAG_MESON_MASS:.8g}", f"{frag_param:.8g}", "9"]
    lines.extend(f"{pt:.10g} {y:.8g}" for pt, y in points)
    lines.append("-1 0")
    return "\n".join(lines) + "\n"


def parse_frag(path: Path) -> dict[tuple[float, float], float]:
    values: dict[tuple[float, float], float] = {}
    for line in path.read_text().splitlines()[1:]:
        cols = line.split()
        if len(cols) < 3:
            continue
        values[point_key(float(cols[0].replace("D", "E")), float(cols[1].replace("D", "E")))] = \
            float(cols[2].replace("D", "E"))
    return values


def trapz2(values: dict[tuple[float, float], float]) -> float:
    total = 0.0
    for ipt, pt in enumerate(PT_VALUES):
        wpt = 0.5 if ipt in (0, len(PT_VALUES) - 1) else 1.0
        for iy, y in enumerate(Y_VALUES):
            wy = 0.5 if iy in (0, len(Y_VALUES) - 1) else 1.0
            total += wpt * wy * values[point_key(pt, y)]
    return total * (PT_VALUES[1] - PT_VALUES[0]) * (Y_VALUES[1] - Y_VALUES[0])


def write_final(out_path: Path, pdf_name: str, pdf_id: int, quark: str, quark_config: dict[str, object],
                values: dict[tuple[float, float], float], grid_file: Path, extra_headers=None,
                mass: float | None = None, ffact: float = 1.0, fren: float = 1.0) -> dict[str, object]:
    """The three-column meson grid with its header; pT = 0 rows are set to 0."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mass = float(quark_config["mass"]) if mass is None else float(mass)
    for y in Y_VALUES:
        values[point_key(0.0, y)] = 0.0
    rows = []
    for pt in PT_VALUES:
        for y in Y_VALUES:
            key = point_key(pt, y)
            if key not in values:
                raise RuntimeError(f"missing fragmented point pT={pt}, y={y}")
            val = values[key]
            if not math.isfinite(val):
                raise RuntimeError(f"non-finite value at pT={pt}, y={y}: {val}")
            rows.append((pt, y, val))
    integral = trapz2(values)
    with out_path.open("w") as out:
        out.write("# FONLL heavy-flavor meson grid\n")
        out.write("# columns: pT y dsigma/dpT/dy\n")
        out.write("# units: GeV 1 pb/GeV\n")
        out.write("# collision: pp\n")
        out.write("# ebeam1_GeV: 7000\n")
        out.write("# ebeam2_GeV: 7000\n")
        out.write("# sqrt_s_GeV: 14000\n")
        out.write(f"# quark: {quark}\n")
        out.write("# final_state: meson\n")
        out.write(f"# final_state_detail: {quark_config['final_state_detail']}\n")
        out.write(f"# heavy_quark_mass_GeV: {mass:.8g}\n")
        out.write(f"# scale: mu0=sqrt(m^2+pT^2); ffact(muF)={ffact:.8g}, fren(muR)={fren:.8g}\n")
        out.write(f"# pdf: {pdf_name}\n")
        out.write(f"# lhapdf_id: {pdf_id}\n")
        out.write("# perturbative_order: FONLL\n")
        out.write(f"# fragmentation_function: {quark_config['fragmentation']}\n")
        out.write(f"# fragmentation_mode: {quark_config['frag_mode']}\n")
        out.write(f"# fragmentation_{quark_config['frag_param_name']}: {quark_config['frag_param']:.8g}\n")
        for key, value in (extra_headers or {}).items():
            out.write(f"# {key}: {value}\n")
        out.write("# fragmentation_frame: y=0\n")
        out.write("# fragmentation_fraction: 1\n")
        out.write("# meson_mass: heavy quark mass\n")
        out.write(f"# internal_quark_grid: {Path(grid_file).name}\n")
        out.write(f"# output_pT_values_GeV: 0..50 step {PT_VALUES[1] - PT_VALUES[0]:.8g}\n")
        out.write(f"# output_y_values: -3..3 step {Y_VALUES[1] - Y_VALUES[0]:.8g}\n")
        out.write("# note_pT0: pT=0 rows are set to 0; fragmfonll evaluates differential points for pT>0\n")
        out.write(f"# trapezoid_integral_pb_y-3to3_pt0to50: {integral:.12e}\n")
        out.write("# pT y dsigma/dpT/dy\n")
        for pt, y, val in rows:
            out.write(f"{pt:.8g} {y:.8g} {val:.12e}\n")
    vals = [r[2] for r in rows]
    return {"path": str(out_path), "rows": len(rows), "pt_min": min(r[0] for r in rows),
            "pt_max": max(r[0] for r in rows), "y_min": min(r[1] for r in rows), "y_max": max(r[1] for r in rows),
            "min_dsigma_dpT_dy_pb_per_GeV": min(vals), "max_dsigma_dpT_dy_pb_per_GeV": max(vals),
            "trapezoid_integral_pb": integral}


def run_fragmentation(ws: Workspace, pdf_key: str, quark: str, run_dir: Path, grid_file: Path, frag_mode: int,
                      frag_param: float, points: list[tuple[float, float]], label: str
                      ) -> tuple[dict[tuple[float, float], float], Path]:
    frag_path = run_dir / "fragmfonll.dat"
    if frag_path.exists():
        frag_path.unlink()
    run_command([str(ws.executable("fragmfonll"))], run_dir,
                frag_input_points(grid_file.name, frag_mode, frag_param, points),
                ws.logs / f"fragmfonll_{pdf_key}_{quark}_{label}.log", ws.env())
    if not frag_path.exists():
        raise RuntimeError(f"expected fragmentation file not found: {frag_path}")
    saved_path = run_dir / f"fragmfonll_{label}.dat"
    shutil.copy2(frag_path, saved_path)
    return parse_frag(frag_path), saved_path


def ensure_charm_calibration_grid(ws: Workspace) -> Path:
    """The CTEQ6.6 charm quark grid at y in {-0.5, 0, 0.5} the feeddown
    weight is fitted on."""
    run_dir = ws.run / "validate_cteq66_charm_y0"
    run_dir.mkdir(parents=True, exist_ok=True)
    grid_file = run_dir / "val_c_.out"
    if grid_file.exists():
        return grid_file
    run_command([str(ws.executable("fonllgridlha"))], run_dir,
                grid_input("val_c_", 10550, 1.5, y_values=[-0.5, 0.0, 0.5]),
                ws.logs / "validate_cteq66_charm_y0_fonllgrid.log", ws.env())
    if not grid_file.exists():
        raise RuntimeError(f"expected calibration grid file not found: {grid_file}")
    return grid_file


def combine_charm_d0_feeddown(direct_values, vector_values, vector_to_direct_weight: float,
                              points: Iterable[tuple[float, float]]) -> dict[tuple[float, float], float]:
    """D0 = (direct + w * sum_c BR_c * D*(pT/z_c)/z_c) / (1 + w * sum_c BR_c)."""
    channels = charm_d0_feeddown_channels()
    branch_sum = sum(float(channel["branching_fraction"]) for channel in channels)
    combined: dict[tuple[float, float], float] = {}
    for pt, y in points:
        feeddown = 0.0
        for channel in channels:
            br = float(channel["branching_fraction"])
            z = float(channel["z_collinear"])
            feeddown += br * vector_values[point_key(pt / z, y)] / z
        direct = direct_values[point_key(pt, y)]
        combined[point_key(pt, y)] = (direct + vector_to_direct_weight * feeddown) / (1.0 + vector_to_direct_weight * branch_sum)
    return combined


def fit_charm_d0_feeddown_weight(ws: Workspace) -> dict[str, object]:
    """Golden-section fit of the D* -> D0 weight to the public CTEQ6.6 D0 points."""
    grid_file = ensure_charm_calibration_grid(ws)
    run_dir = grid_file.parent
    reference_points = [(pt, 0.0) for pt in sorted(PUBLIC_CHARM_D0_CTEQ66_Y0)]
    parent_points = charm_feeddown_parent_points(reference_points)
    direct_values, _ = run_fragmentation(ws, "cteq66", "charm", run_dir, grid_file, 5, 0.1, reference_points,
                                         "calibration_direct_pseudoscalar")
    vector_values, _ = run_fragmentation(ws, "cteq66", "charm", run_dir, grid_file, 4, 0.1, parent_points,
                                         "calibration_dstar_vector")

    def objective(weight: float) -> float:
        combined = combine_charm_d0_feeddown(direct_values, vector_values, weight, reference_points)
        return sum(math.log(combined[point_key(pt, 0.0)] / public) ** 2
                   for pt, public in PUBLIC_CHARM_D0_CTEQ66_Y0.items())

    lo, hi = 0.0, 5.0
    for _ in range(100):
        left = lo + (hi - lo) / 3.0
        right = hi - (hi - lo) / 3.0
        if objective(left) < objective(right):
            hi = right
        else:
            lo = left
    weight = 0.5 * (lo + hi)
    combined = combine_charm_d0_feeddown(direct_values, vector_values, weight, reference_points)
    point_checks, max_abs_rel = [], 0.0
    for pt, public in sorted(PUBLIC_CHARM_D0_CTEQ66_Y0.items()):
        local = combined[point_key(pt, 0.0)]
        rel = (local - public) / public
        max_abs_rel = max(max_abs_rel, abs(rel))
        point_checks.append({"pt": pt, "y": 0.0, "local": local, "public": public, "relative_difference": rel})
    return {"grid": str(grid_file), "vector_to_direct_weight": weight,
            "max_abs_relative_difference": max_abs_rel, "points": point_checks}


def load_or_fit_charm_feeddown_weight(ws: Workspace) -> dict[str, object]:
    """Calibrate the feeddown weight once and cache it.

    The weight is anchored to fixed public CTEQ6.6 references, so it is a
    property of the fragmentation model, not of the scale/PDF/mass
    variation; every charm variation reuses this single calibration.
    """
    cache_path = ws.out / "charm_feeddown_calibration.json"
    if cache_path.exists():
        return json.loads(cache_path.read_text())
    calibration = fit_charm_d0_feeddown_weight(ws)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(calibration, indent=2) + "\n")
    return calibration


def short_pdf(pdf_key: str) -> str:
    pdf_entry = PDFS[pdf_key]
    if len(pdf_entry) == 3:
        return pdf_entry[2]
    return "nnpdf40_nlo_as_01180" if pdf_key == "nlo" else "nnpdf40_nnlo_as_01180"


def grid_filename(pdf_key: str, tag: str, quark: str) -> str:
    return f"fonll_pp14tev_{short_pdf(pdf_key)}_fonll_meson_dsdpTdy_pt0-50_y-3to3_{tag}_{quark}.dat"


def generate(ws: Workspace, pdf_key: str, quark: str, variation: Variation, reuse_existing_grids: bool,
             grid_workers: int, feeddown_calibration=None) -> dict[str, object]:
    """One meson grid: quark grid (unless reused), fragmentation, header, summary."""
    pdf_name, base_pdf_id = PDFS[pdf_key][:2]
    quark_config = QUARKS[quark]
    mass = variation.mass
    lhaid = int(base_pdf_id) + int(variation.pdf_member)
    ffact, fren = variation.muF, variation.muR
    out_path = ws.out / grid_filename(pdf_key, variation.tag, quark)
    run_dir = ws.run / f"{pdf_key}_{quark}_{variation.tag}"
    run_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"{pdf_key[:2]}_{quark[:1]}_"
    grid_file = run_dir / f"{prefix}.out"
    stale = [run_dir / f"{prefix}.outlog", run_dir / f"{prefix}fonll.log", run_dir / f"{prefix}hdml.tmp",
             run_dir / f"{prefix}hdmv.tmp", run_dir / f"{prefix}hdrs.tmp", run_dir / "fragmfonll.dat",
             run_dir / "fragfonll.log"]
    if not reuse_existing_grids:
        stale.append(grid_file)
    for path in stale:
        if path.exists():
            path.unlink()

    label = f"[{pdf_key} {quark} {variation.tag}]"
    if reuse_existing_grids and grid_file.exists():
        print(f"{label} reusing existing quark grid {grid_file}")
    else:
        print(f"{label} building quark grid")
        build_quark_grid(ws, run_dir, grid_file, prefix, lhaid, mass, pdf_key, quark, grid_workers,
                         ffact=ffact, fren=fren, log_label=variation.tag)
    if not grid_file.exists():
        raise RuntimeError(f"expected grid file not found: {grid_file}")

    points = observed_points()
    extra_headers: dict[str, object] = {"lhapdf_member": variation.pdf_member, "variation_kind": variation.kind,
                                        "variation_tag": variation.tag}
    weight = None
    if quark_config.get("feeddown") == "public_d0":
        if feeddown_calibration is None:
            feeddown_calibration = load_or_fit_charm_feeddown_weight(ws)
        weight = float(feeddown_calibration["vector_to_direct_weight"])
        print(f"{label} applying direct pseudoscalar charm fragmentation")
        direct_values, _ = run_fragmentation(ws, pdf_key, f"{quark}_{variation.tag}", run_dir, grid_file,
                                             int(quark_config["frag_mode"]), float(quark_config["frag_param"]),
                                             points, "direct_pseudoscalar")
        print(f"{label} applying D* vector fragmentation for feeddown")
        vector_values, _ = run_fragmentation(ws, pdf_key, f"{quark}_{variation.tag}", run_dir, grid_file, 4,
                                             float(quark_config["frag_param"]),
                                             charm_feeddown_parent_points(points), "dstar_vector_feeddown")
        values = combine_charm_d0_feeddown(direct_values, vector_values, weight, points)
        extra_headers.update({
            "charm_feeddown_model": "collinear two-body D* -> D0 feeddown, calibrated to public FONLL CTEQ6.6 D0",
            "charm_feeddown_vector_to_direct_weight": f"{weight:.12g}",
            "charm_feeddown_branching_sum": f"{sum(float(c['branching_fraction']) for c in charm_d0_feeddown_channels()):.12g}",
            "charm_feeddown_cteq66_reference": "public FONLL v1.3.2 CTEQ6.6, queried 2026-05-29",
            "charm_feeddown_cteq66_fit_max_abs_rel_diff": f"{float(feeddown_calibration['max_abs_relative_difference']):.6e}",
        })
    else:
        print(f"{label} applying meson fragmentation")
        values, _ = run_fragmentation(ws, pdf_key, f"{quark}_{variation.tag}", run_dir, grid_file,
                                      int(quark_config["frag_mode"]), float(quark_config["frag_param"]), points,
                                      "central")

    summary = write_final(out_path, pdf_name, lhaid, quark, quark_config, values, grid_file,
                          extra_headers=extra_headers, mass=mass, ffact=ffact, fren=fren)
    summary.update({"pdf_key": pdf_key, "quark": quark, "variation_kind": variation.kind,
                    "variation_tag": variation.tag, "muR": variation.muR, "muF": variation.muF, "lhapdf_id": lhaid,
                    "lhapdf_member": variation.pdf_member, "heavy_quark_mass_GeV": mass,
                    "raw_grid_reused": reuse_existing_grids})
    if weight is not None:
        summary["charm_feeddown_vector_to_direct_weight"] = weight
    if not reuse_existing_grids:
        summary["raw_grid_workers"] = grid_workers
    return summary


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fonll_version(ws: Workspace) -> str:
    readme = ws.fonll / "README"
    if readme.exists():
        for line in readme.read_text(errors="ignore").splitlines():
            stripped = line.strip()
            if "version" in stripped.lower() or stripped.lower().startswith("fonll"):
                return stripped
    return "FONLL (version unrecorded; see the README of the FONLL tree)"


PATCH_NOTES = ["misc1/fragmfonll.f: BCFY vector/pseudoscalar modes 4 and 5/8",
               "misc1/fonllgrid.f: rapidity capacity raised to 120 nodes"]


def manifest_document(ws: Workspace, entries: list[dict], pdf_name: str, base_pdf_id: int, pdf_members: list[int],
                      scale_points, mass_variations, feeddown_calibration, note: str | None = None) -> dict:
    """The shared shape of ``variation_manifest.json`` (campaign and snapshot)."""
    manifest = {"generated_unix": time.time(), "generated_iso": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    if note:
        manifest["note"] = note
    manifest.update({
        "fonll_version": fonll_version(ws), "patches": PATCH_NOTES,
        "collision": "pp", "ebeam1_GeV": 7000, "ebeam2_GeV": 7000, "sqrt_s_GeV": 14000,
        "pdf_set": pdf_name, "pdf_base_lhapdf_id": base_pdf_id, "pdf_members": list(pdf_members),
        "scale_points_muR_muF": scale_points, "mass_variations_GeV": mass_variations,
        "grid_bounds": {"pt_min_GeV": PT_VALUES[0], "pt_max_GeV": PT_VALUES[-1], "y_min": Y_VALUES[0],
                        "y_max": Y_VALUES[-1], "n_pt": len(PT_VALUES), "n_y": len(Y_VALUES)},
        "charm_feeddown": feeddown_calibration, "grids": entries,
    })
    return manifest


def write_manifest(ws: Workspace, summaries: list[dict[str, object]], pdf_key: str, pdf_members: list[int],
                   include_scale: bool, include_mass: bool, feeddown_calibration) -> Path:
    pdf_name, base_pdf_id = PDFS[pdf_key][:2]
    entries = []
    for summary in sorted(summaries, key=lambda s: (s["quark"], s["variation_tag"])):
        path = Path(str(summary["path"]))
        entries.append({key: summary[key] for key in ("quark", "variation_kind", "variation_tag", "muR", "muF",
                                                       "lhapdf_id", "lhapdf_member", "heavy_quark_mass_GeV")}
                       | {"path": str(path), "sha256": sha256_file(path), "rows": summary["rows"],
                          "trapezoid_integral_pb": summary["trapezoid_integral_pb"]})
    manifest = manifest_document(ws, entries, pdf_name, base_pdf_id, pdf_members,
                                 SCALE_POINTS_7 if include_scale else [(1.0, 1.0)],
                                 MASS_VARIATIONS if include_mass else {}, feeddown_calibration)
    manifest_path = ws.out / "variation_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest_path


def parse_pdf_members(spec: str) -> list[int]:
    """'0-100', '0,1,5' ... ; member 0 is always included."""
    if spec.strip().lower() in {"", "none"}:
        return [0]
    members: list[int] = []
    for part in spec.split(","):
        part = part.strip()
        if "-" in part:
            lo, hi = part.split("-", 1)
            members.extend(range(int(lo), int(hi) + 1))
        else:
            members.append(int(part))
    if 0 not in members:
        members.insert(0, 0)
    return sorted(dict.fromkeys(members))


def run_campaign(ws: Workspace, pdf_key: str, quarks: list[str], pdf_members: list[int], include_scale: bool,
                 include_mass: bool, max_parallel: int, grid_workers: int, reuse_existing_grids: bool) -> Path:
    feeddown_calibration = load_or_fit_charm_feeddown_weight(ws) if "charm" in quarks else None
    tasks = [(pdf_key, quark, variation) for quark in quarks
             for variation in enumerate_variations(quark, pdf_members, include_scale, include_mass)]
    print(f"campaign: {len(tasks)} variation grids ({', '.join(quarks)}; scale={include_scale}; "
          f"pdf_members={len(pdf_members)}; mass={include_mass}); max_parallel={max_parallel}, "
          f"grid_workers={grid_workers}")
    summaries: list[dict[str, object]] = []
    failures = []

    def one(task):
        pk, q, var = task
        return generate(ws, pk, q, var, reuse_existing_grids, grid_workers,
                        feeddown_calibration=feeddown_calibration if q == "charm" else None)

    if max_parallel <= 1:
        for task in tasks:
            summaries.append(one(task))
    else:
        with ThreadPoolExecutor(max_workers=max_parallel) as executor:
            futures = {executor.submit(one, task): task for task in tasks}
            for future in as_completed(futures):
                pk, q, var = futures[future]
                try:
                    summaries.append(future.result())
                except Exception as exc:  # noqa: BLE001 - report and continue
                    failures.append(((pk, q, var), str(exc)))
                    print(f"[FAILED] {pk} {q} {var.tag}: {exc}")
    manifest_path = write_manifest(ws, summaries, pdf_key, pdf_members, include_scale, include_mass,
                                   feeddown_calibration)
    print(f"wrote {manifest_path} ({len(summaries)} grids)")
    if failures:
        print(f"WARNING: {len(failures)} variation grids failed:")
        for (pk, q, var), msg in failures:
            print(f"  - {pk} {q} {var.tag}: {msg}")
    return manifest_path


def add_workspace_options(ap: argparse.ArgumentParser) -> None:
    ap.add_argument("--fonll", type=Path, default=None, help="FONLL tree with Linux/{fonllgridlha,fragmfonll} "
                                                            "(default: GRENDEL_FONLL_DIR or third_party/fonll)")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="grid output directory (default: GRENDEL_FONLL_GRID_DIR or <work>/fonll_grids/output)")
    ap.add_argument("--run-dir", type=Path, default=None, help="scratch for the FONLL runs (default: beside --out-dir)")
    ap.add_argument("--log-dir", type=Path, default=None, help="per-command logs (default: beside --out-dir)")


def workspace_from(args) -> Workspace:
    return Workspace.resolve(args.fonll, args.out_dir, args.run_dir, args.log_dir)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_workspace_options(ap)
    ap.add_argument("--pdf", action="append", choices=sorted(PDFS),
                    help="PDF key; may be repeated (default: nlo, plus nnlo with --include-nnlo)")
    ap.add_argument("--include-nnlo", action="store_true", help="also generate the NNLO comparison grids")
    ap.add_argument("--quark", action="append", choices=sorted(QUARKS),
                    help="quark flavour; may be repeated (default: bottom and charm)")
    ap.add_argument("--reuse-existing-grids", action="store_true",
                    help="skip fonllgridlha when the raw quark grid already exists")
    ap.add_argument("--grid-workers", type=int, default=1, help="parallel rapidity chunks per quark grid")
    ap.add_argument("--campaign", action="store_true",
                    help="run the variation campaign (scale + PDF members + mass) and write variation_manifest.json")
    ap.add_argument("--scale-variations", action="store_true", help="include the 7-point scale variations")
    ap.add_argument("--pdf-members", default="0", help="PDF members, e.g. '0-100' or '0,1,5' (default: central only)")
    ap.add_argument("--mass-variations", action="store_true", help="include the heavy-quark-mass variations")
    ap.add_argument("--max-parallel", type=int, default=1, help="variation grids computed concurrently")
    ap.add_argument("--compress-logs", action="store_true", help="gzip each per-command log on success")
    args = ap.parse_args(argv)
    if args.grid_workers < 1:
        ap.error("--grid-workers must be at least 1")
    if args.max_parallel < 1:
        ap.error("--max-parallel must be at least 1")
    global _COMPRESS_LOGS
    _COMPRESS_LOGS = args.compress_logs

    ws = workspace_from(args)
    ws.out.mkdir(parents=True, exist_ok=True)
    requested = list(dict.fromkeys((args.pdf or ["nlo"]) + (["nnlo"] if args.include_nnlo else [])))
    quarks = args.quark or ["bottom", "charm"]

    if args.campaign:
        narrowed = args.scale_variations or args.mass_variations or args.pdf_members != "0"
        pdf_members = parse_pdf_members(args.pdf_members)
        if narrowed:
            include_scale, include_mass = args.scale_variations, args.mass_variations
        else:
            # Bare --campaign: the full default set (scale + all 100 members + mass).
            include_scale = include_mass = True
            pdf_members = list(range(0, 101))
        for pdf_key in requested:
            run_campaign(ws, pdf_key, quarks, pdf_members, include_scale, include_mass, args.max_parallel,
                         args.grid_workers, args.reuse_existing_grids)
        return 0

    feeddown_calibration = load_or_fit_charm_feeddown_weight(ws) if "charm" in quarks else None
    summaries = [generate(ws, pdf_key, quark, central_variation(quark), args.reuse_existing_grids, args.grid_workers,
                          feeddown_calibration=feeddown_calibration if quark == "charm" else None)
                 for pdf_key in requested for quark in quarks]
    summary_name = ("generation_summary.json" if requested == ["nlo"] and quarks == ["bottom", "charm"]
                    else f"generation_summary_{'_'.join(requested + quarks)}.json")
    summary_path = ws.out / summary_name
    summary_path.write_text(json.dumps(summaries, indent=2) + "\n")
    print(f"wrote {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
