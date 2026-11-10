"""Combine the variation grids into uncertainty envelopes.

    python -m grendel.production.fonll_grids.combine variations [--out-dir DIR]
    python -m grendel.production.fonll_grids.combine alphas     [--out-dir DIR]

``variations`` reads ``variation_manifest.json``, partitions each quark's
grids by axis and writes envelope grids under ``<out-dir>/envelopes`` in
the same three-column format:

* scale: pointwise max/min over the 7-point (muR, muF) set;
* pdf:   NNPDF Monte-Carlo prescription, central +/- std (ddof=1) over the
         replica members (member 0 is the replica mean, not a replica, and
         is excluded); the replica mean is also written;
* mass:  pointwise max/min over the m_b / m_c up/down set;
* combined: central +/- the three deviations in uncorrelated quadrature.

``alphas`` writes the PDF4LHC strong-coupling band, half the difference of
the ``NNPDF40_nlo_as_01170`` / ``_01190`` companion central grids applied
symmetrically around the ``as_01180`` central, to the same directory with
its own ``alphas_manifest.json``. The combined envelope does NOT include
alpha_s; add the two in quadrature for the total.

Envelopes are pointwise constructions, not physical cross sections. They
carry an ``envelope_band`` header and the production sampler refuses them
(``grendel.production.fonll.fonll_parser``); only coherent individual grids
(one scale point, one replica, one mass) may be sampled.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

from . import generate as gen

PT_VALUES = gen.PT_VALUES
Y_VALUES = gen.Y_VALUES
QUARKS = gen.QUARKS
ATOL = 1e-9
ALPHAS_SETS = {"down": ("nlo_as_01170", 0.117), "central": ("nlo", 0.118), "up": ("nlo_as_01190", 0.119)}


def _close(a: float, b: float) -> bool:
    return abs(float(a) - float(b)) < ATOL


def reference_grid_coords() -> tuple[np.ndarray, np.ndarray]:
    """Canonical (pT, y) column order: pT outermost, y innermost, ascending."""
    pt = np.repeat(np.array(PT_VALUES, dtype=float), len(Y_VALUES))
    y = np.tile(np.array(Y_VALUES, dtype=float), len(PT_VALUES))
    return pt, y


def load_grid_column(path: Path, ref_pt: np.ndarray, ref_y: np.ndarray) -> np.ndarray:
    """The dsigma column of a grid in the canonical ordering."""
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")
    data = np.loadtxt(path, comments="#")
    if data.shape != (len(ref_pt), 3):
        raise ValueError(f"{path}: expected {len(ref_pt)}x3 grid, got {data.shape}")
    if not (np.allclose(data[:, 0], ref_pt, atol=1e-4) and np.allclose(data[:, 1], ref_y, atol=1e-4)):
        raise ValueError(f"{path}: grid coordinates do not match the canonical pT/y ordering")
    col = data[:, 2]
    if not np.all(np.isfinite(col)):
        raise ValueError(f"{path}: non-finite dsigma values")
    return col


def trapz2(column: np.ndarray) -> float:
    grid = column.reshape(len(PT_VALUES), len(Y_VALUES))
    return float(np.trapezoid(np.trapezoid(grid, Y_VALUES, axis=1), PT_VALUES))


def partition_entries(entries: list[dict], quark: str) -> dict[str, object]:
    """central / scale / pdf / mass subsets of one quark's manifest entries."""
    default_mass = float(QUARKS[quark]["mass"])
    q = [e for e in entries if e["quark"] == quark]

    def is_central_scale(e):
        return _close(e["muR"], 1.0) and _close(e["muF"], 1.0)

    def is_central_mass(e):
        return _close(e["heavy_quark_mass_GeV"], default_mass)

    central = next((e for e in q if e["lhapdf_member"] == 0 and is_central_scale(e) and is_central_mass(e)), None)
    if central is None:
        raise ValueError(f"{quark}: no central grid (member 0, scale 1,1, mass {default_mass}) in manifest")
    return {"central": central,
            "scale": sorted((e for e in q if e["lhapdf_member"] == 0 and is_central_mass(e)),
                            key=lambda e: (e["muR"], e["muF"])),
            "pdf": sorted((e for e in q if is_central_scale(e) and is_central_mass(e)),
                          key=lambda e: e["lhapdf_member"]),
            "mass": sorted((e for e in q if e["lhapdf_member"] == 0 and is_central_scale(e)),
                           key=lambda e: e["heavy_quark_mass_GeV"]),
            "default_mass": default_mass}


def write_envelope_grid(out_path: Path, quark: str, band: str, column: np.ndarray, ref_pt: np.ndarray,
                        ref_y: np.ndarray, provenance: dict, title: str = "variation envelope") -> dict:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    column = np.maximum(column, 0.0)
    integral = trapz2(column)
    with out_path.open("w") as out:
        out.write(f"# FONLL heavy-flavor meson grid ({title})\n")
        out.write("# columns: pT y dsigma/dpT/dy\n")
        out.write("# units: GeV 1 pb/GeV\n")
        out.write("# collision: pp\n")
        out.write("# sqrt_s_GeV: 14000\n")
        out.write(f"# quark: {quark}\n")
        out.write(f"# envelope_band: {band}\n")
        for key, value in provenance.items():
            out.write(f"# {key}: {value}\n")
        out.write(f"# trapezoid_integral_pb_y-3to3_pt0to50: {integral:.12e}\n")
        out.write("# pT y dsigma/dpT/dy\n")
        for pt, y, val in zip(ref_pt, ref_y, column):
            out.write(f"{pt:.8g} {y:.8g} {val:.12e}\n")
    return {"band": band, "path": str(out_path), "sha256": gen.sha256_file(out_path),
            "trapezoid_integral_pb": integral}


def combine_quark(quark: str, manifest: dict, env_dir: Path) -> list[dict]:
    parts = partition_entries(manifest["grids"], quark)
    ref_pt, ref_y = reference_grid_coords()
    central_col = load_grid_column(Path(parts["central"]["path"]), ref_pt, ref_y)
    stem = (f"fonll_pp14tev_{str(manifest.get('pdf_set', '')).lower().replace('.', '')}_fonll_meson_dsdpTdy_"
            f"pt0-50_y-3to3_envelope")
    written: list[dict] = []

    def write(band_tag: str, band: str, column, prov):
        written.append(write_envelope_grid(env_dir / f"{stem}_{band_tag}_{quark}.dat", quark, band, column,
                                           ref_pt, ref_y, prov))

    write("central", "central", central_col, {"source": parts["central"]["path"]})

    scaleup = scaledn = None
    if len(parts["scale"]) >= 2:
        cols = np.stack([load_grid_column(Path(e["path"]), ref_pt, ref_y) for e in parts["scale"]])
        scaleup, scaledn = cols.max(axis=0), cols.min(axis=0)
        prov = {"n_grids": len(parts["scale"]),
                "scale_points_muR_muF": "; ".join(f"({e['muR']:g},{e['muF']:g})" for e in parts["scale"])}
        write("scaleup", "scale_up", scaleup, prov)
        write("scaledn", "scale_dn", scaledn, prov)
    else:
        print(f"[{quark}] scale band skipped (only {len(parts['scale'])} grid)")

    pdf_sigma = None
    replicas = [e for e in parts["pdf"] if int(e["lhapdf_member"]) >= 1]
    if len(replicas) >= 2:
        cols = np.stack([load_grid_column(Path(e["path"]), ref_pt, ref_y) for e in replicas])
        pdf_mean = cols.mean(axis=0)
        pdf_sigma = cols.std(axis=0, ddof=1)
        prov = {"n_replicas": len(replicas),
                "prescription": "NNPDF Monte-Carlo: central +/- std(ddof=1) over replica members >= 1 "
                                "(member 0 excluded); replica mean also written"}
        write("pdfmean", "pdf_mean", pdf_mean, prov)
        write("pdfup", "pdf_up", central_col + pdf_sigma, prov)
        write("pdfdn", "pdf_dn", central_col - pdf_sigma, prov)
    else:
        print(f"[{quark}] pdf band skipped (only {len(replicas)} replica members)")

    massup = massdn = None
    if len(parts["mass"]) >= 2:
        cols = np.stack([load_grid_column(Path(e["path"]), ref_pt, ref_y) for e in parts["mass"]])
        massup, massdn = cols.max(axis=0), cols.min(axis=0)
        prov = {"n_grids": len(parts["mass"]),
                "masses_GeV": "; ".join(f"{e['heavy_quark_mass_GeV']:g}" for e in parts["mass"])}
        write("massup", "mass_up", massup, prov)
        write("massdn", "mass_dn", massdn, prov)
    else:
        print(f"[{quark}] mass band skipped (only {len(parts['mass'])} grid)")

    components = []
    up_var = np.zeros_like(central_col)
    dn_var = np.zeros_like(central_col)
    if scaleup is not None:
        up_var += (scaleup - central_col) ** 2
        dn_var += (central_col - scaledn) ** 2
        components.append("scale")
    if pdf_sigma is not None:
        up_var += pdf_sigma ** 2
        dn_var += pdf_sigma ** 2
        components.append("pdf")
    if massup is not None:
        up_var += (massup - central_col) ** 2
        dn_var += (central_col - massdn) ** 2
        components.append("mass")
    if components:
        prov = {"components": "+".join(components),
                "treatment": "uncorrelated quadrature of per-axis deviations from central"}
        write("combup", "combined_up", central_col + np.sqrt(up_var), prov)
        write("combdn", "combined_dn", central_col - np.sqrt(dn_var), prov)
    return written


def combine_variations(ws: gen.Workspace, quarks: list[str]) -> Path:
    manifest_path = ws.out / "variation_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"{manifest_path} not found; run generate --campaign (or snapshot) first")
    manifest = json.loads(manifest_path.read_text())
    available = {e["quark"] for e in manifest["grids"]}
    env_dir = ws.out / "envelopes"
    env_dir.mkdir(parents=True, exist_ok=True)
    envelopes = {}
    for quark in quarks:
        if quark not in available:
            print(f"[{quark}] no grids in manifest; skipping")
            continue
        envelopes[quark] = combine_quark(quark, manifest, env_dir)
    summary_path = env_dir / "envelope_manifest.json"
    summary_path.write_text(json.dumps({"generated_iso": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                                        "source_manifest": str(manifest_path), "envelopes": envelopes},
                                       indent=2) + "\n")
    print(f"wrote {summary_path}")
    for quark, grids in envelopes.items():
        print(f"  {quark}: {len(grids)} envelope grids -> {', '.join(g['band'] for g in grids)}")
    return summary_path


def combine_alphas_quark(ws: gen.Workspace, quark: str, env_dir: Path) -> dict:
    ref_pt, ref_y = reference_grid_coords()
    cols = {}
    for name, (pdf_key, _) in ALPHAS_SETS.items():
        path = ws.out / gen.grid_filename(pdf_key, "central", quark)
        if not path.exists():
            raise FileNotFoundError(f"{path} not found; generate the alpha_s companion grids first "
                                    f"(generate --pdf nlo_as_01170 --pdf nlo_as_01190 --quark {quark})")
        cols[name] = load_grid_column(path, ref_pt, ref_y)
    halfdiff = 0.5 * (cols["up"] - cols["down"])     # signed: as = 0.119 minus as = 0.117
    band = np.abs(halfdiff)
    stem = f"fonll_pp14tev_{gen.short_pdf('nlo')}_fonll_meson_dsdpTdy_pt0-50_y-3to3_alphas"
    prov = {"prescription": "PDF4LHC alpha_s: 0.5*(sigma(as=0.119) - sigma(as=0.117)) around as=0.118 central",
            "alphas_down": ALPHAS_SETS["down"][1], "alphas_central": ALPHAS_SETS["central"][1],
            "alphas_up": ALPHAS_SETS["up"][1]}
    written = [write_envelope_grid(env_dir / f"{stem}up_{quark}.dat", quark, "alphas_up", cols["central"] + band,
                                   ref_pt, ref_y, prov, title="alpha_s envelope"),
               write_envelope_grid(env_dir / f"{stem}dn_{quark}.dat", quark, "alphas_dn", cols["central"] - band,
                                   ref_pt, ref_y, prov, title="alpha_s envelope"),
               write_envelope_grid(env_dir / f"{stem}halfdiff_{quark}.dat", quark, "alphas_halfdiff", band,
                                   ref_pt, ref_y, prov, title="alpha_s envelope")]
    int_c, int_band = trapz2(cols["central"]), trapz2(band)
    return {"quark": quark, "alphas_down": ALPHAS_SETS["down"][1], "alphas_central": ALPHAS_SETS["central"][1],
            "alphas_up": ALPHAS_SETS["up"][1], "central_integral_pb": int_c, "halfdiff_integral_pb": int_band,
            "integrated_alphas_uncertainty_fraction": (int_band / int_c) if int_c else float("nan"),
            "envelopes": written}


def combine_alphas(ws: gen.Workspace, quarks: list[str]) -> Path:
    env_dir = ws.out / "envelopes"
    env_dir.mkdir(parents=True, exist_ok=True)
    per_quark = {}
    for quark in quarks:
        result = combine_alphas_quark(ws, quark, env_dir)
        per_quark[quark] = result
        print(f"[{quark}] alpha_s band: integrated +/-{result['integrated_alphas_uncertainty_fraction'] * 100:.2f}% "
              f"of central; {len(result['envelopes'])} envelope grids")
    summary_path = env_dir / "alphas_manifest.json"
    summary_path.write_text(json.dumps({
        "generated_iso": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "prescription": "PDF4LHC alpha_s half-difference of NNPDF40_nlo_as_01170/01190 around as_01180",
        "quarks": per_quark}, indent=2) + "\n")
    print(f"wrote {summary_path}")
    return summary_path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=("variations", "alphas"))
    gen.add_workspace_options(ap)
    ap.add_argument("--quark", action="append", choices=sorted(QUARKS),
                    help="quark to combine; may be repeated (default: bottom and charm)")
    args = ap.parse_args(argv)
    ws = gen.workspace_from(args)
    quarks = args.quark or ["bottom", "charm"]
    if args.what == "variations":
        combine_variations(ws, quarks)
    else:
        combine_alphas(ws, quarks)
    return 0


if __name__ == "__main__":
    sys.exit(main())
