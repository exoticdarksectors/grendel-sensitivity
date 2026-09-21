#!/usr/bin/env python3
"""Verify and export the pinned SensCalc 2501 LHC production inputs.

The production-probability and pregenerated-fragmentation files are Wolfram
MX dumps.  This wrapper pins every source used by the BC10 production audit,
runs the Wolfram decoder in a staging directory, and installs the result only
after its manifest and hashes are complete.

Example:

    python alp_fermion/tools/export_senscalc_2501_production.py /path/to/SensCalc
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from grendel.io.paths import work_dir

try:
    from .export_senscalc_2501 import (
        InputError,
        SENSCALC_COMMIT,
        SENSCALC_TAG,
        git_head,
        sha256,
        source_mx_system_ids,
        wolfram_environment,
    )
except ImportError:  # Direct execution from alp_fermion/tools.
    from export_senscalc_2501 import (  # type: ignore[no-redef]
        InputError,
        SENSCALC_COMMIT,
        SENSCALC_TAG,
        git_head,
        sha256,
        source_mx_system_ids,
        wolfram_environment,
    )


SOURCE_FILES = {
    "sensitivity_notebook": (
        Path("3. ALP-fermion sensitivity.nb"),
        "48883fb7ecce956bd73965709ecf44d0cae286e637067c672c1566d736072ef3",
    ),
    "generic_notebook": (
        Path("codes/generic.nb"),
        "95cc528268e41b273182377275bf9eeb00d8cc1cb92d149ddeaa9e72967f0d4a",
    ),
    "production_notebook": (
        Path("codes/LLP distribution/prod-pheno-ALP-fermion.nb"),
        "e1a7b91a4e82421081d165561a1b7f23e90aefe451efca32a5a3459b1ac0d0f7",
    ),
    "experiments": (
        Path("codes/experiments.nb"),
        "3e60809105f49d8c5274de9a352cd0063bac9cb8feba31188e5e49ca6e468bdf",
    ),
    "coefficients": (
        Path("phenomenology/ALP-fermion/Production probabilities")
        / "Coefficients-ALP-fermion-Lambda=1.-TeV.mx",
        "e6885569ee8d33fb5310723ba2f020bf8b4e3506c2b5ad258b5230f9aa981f43",
    ),
    "fragmentation_probability": (
        Path("phenomenology/ALP-fermion/Production probabilities")
        / "ProductionProbability-Fragmentation-ALP-fermion.mx",
        "4804205cfb14c0d40322fcabb7ec65843023ee44e91f6918603bca0e8835e9eb",
    ),
    "light_meson_decays": (
        Path("phenomenology/ALP-fermion/Production probabilities")
        / "BrRatios-Msquared-LightMesonDecays-ALP-fermion-Lambda=1.-TeV.mx",
        "c7468dcf8863a8ac182c2720d34d12b2dd12693700aec032460cd5d678985f71",
    ),
    "bremsstrahlung_lhc": (
        Path("phenomenology/ALP-fermion/Production probabilities")
        / "ProdProb_ALP-fermion_Bremsstrahlung-AP_LHC.m",
        "435022c25a2bb908ffd5c4b12a5645a0d5aa7f6914c49aab5e33fe0b8b60e3d2",
    ),
    "drell_yan_lhc": (
        Path("phenomenology/ALP-fermion/Production probabilities")
        / "sigmaDrellYan_ALP-fermion_LHC.txt",
        "12e013f6df9a65c07514e75d8f1c47d53b53ebc72b21624ab13eb70afaed185d",
    ),
    "fragmentation_lhc_grid": (
        Path("spectra/New physics particles spectra/ALP-fermion/Pregenerated")
        / "DoubleDistr_ALP-fermion_Fragmentation_LHC.m",
        "43116651d526358c446af54a16d43ca9b2233fb42ff7ab4bae6946fe21477a79",
    ),
    "eta_lhc_spectrum": (
        Path("spectra/SM particles/DoubleDistr_LHC_Eta.txt"),
        "a7e2323c721ce2f6b4b34571b9eaa956ce2cbb5ec8d5c437a2482cafdbfbf39c",
    ),
    "eta_prime_lhc_spectrum": (
        Path("spectra/SM particles/DoubleDistr_LHC_EtaPr.txt"),
        "dfe394e4901f9d3ddf1e2e3881abe2f8fe550a1cbd0bc7a759f87b79fdf8882d",
    ),
    "omega_lhc_spectrum": (
        Path("spectra/SM particles/DoubleDistr_LHC_Omega.txt"),
        "86f8a67260340de8091936a602eed44f33c5d78200afbaca5b5943fb5849eb02",
    ),
    "rho_charged_lhc_spectrum": (
        Path("spectra/SM particles/DoubleDistr_LHC_RhoCh.txt"),
        "b2265c828c8b98762e6b7d99973ed6a97f5f66b01c2cab87708236f9ba4738ff",
    ),
    "ks_lhc_spectrum": (
        Path("spectra/SM particles/DoubleDistr_LHC_KS.txt"),
        "53aafc2e551fed9d01698ba37c7f230e6c5e66c30cc41e73976cd06dd21784ef",
    ),
}

EXPECTED_OUTPUT_FILES = (
    "production_coupling_coefficients_bnt.csv",
    "production_coupling_coefficients_metadata.json",
    "fragmentation_probability_coefficients_bnt.csv",
    "fragmentation_probability_metadata.json",
    "meson_branching_coefficients_bnt.csv",
    "meson_decay_channels.json",
    "meson_decay_matrix_elements.json",
    "fragmentation_lhc_grid.csv",
    "drell_yan_lhc_cross_section_coefficients.csv",
    "drell_yan_lhc_metadata.json",
    "bremsstrahlung_lhc_metadata.json",
    "PRODUCTION_EXPORT_MANIFEST.json",
)

DEFAULT_OUTPUT_DIR = work_dir() / "bc10" / "senscalc_2501_production"
EXPORTER = Path(__file__).with_suffix(".wls")


def verify_sources(repo: Path) -> dict[str, Path]:
    """Return exact, hash-verified production source paths."""
    repo = repo.resolve()
    head = git_head(repo)
    if head != SENSCALC_COMMIT:
        raise InputError(
            f"SensCalc HEAD is {head}, expected {SENSCALC_TAG} "
            f"({SENSCALC_COMMIT})"
        )

    verified = {}
    for key, (relative_path, expected_hash) in SOURCE_FILES.items():
        path = repo / relative_path
        if not path.is_file():
            raise InputError(f"missing SensCalc input: {path}")
        actual_hash = sha256(path)
        if actual_hash != expected_hash:
            raise InputError(
                f"SHA-256 mismatch for {relative_path}: {actual_hash}; "
                f"expected {expected_hash}"
            )
        verified[key] = path
    return verified


def validate_export(export_dir: Path) -> None:
    """Reject incomplete or internally inconsistent Wolfram output."""
    missing = [
        name for name in EXPECTED_OUTPUT_FILES
        if not (export_dir / name).is_file()
    ]
    if missing:
        raise InputError(f"Wolfram export is missing: {', '.join(missing)}")

    manifest_path = export_dir / "PRODUCTION_EXPORT_MANIFEST.json"
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise InputError(
            "PRODUCTION_EXPORT_MANIFEST.json is not valid JSON"
        ) from exc

    if manifest.get("senscalc_tag") != SENSCALC_TAG:
        raise InputError("export manifest has the wrong SensCalc tag")
    if manifest.get("senscalc_commit") != SENSCALC_COMMIT:
        raise InputError("export manifest has the wrong SensCalc commit")
    if manifest.get("phenomenology") != "arXiv:2501.04525":
        raise InputError("export manifest has the wrong phenomenology source")
    if not isinstance(manifest.get("wolfram_system_id"), str):
        raise InputError("export manifest has no Wolfram system ID")
    expected_source_hashes = {
        key: expected_hash
        for key, (_, expected_hash) in SOURCE_FILES.items()
    }
    if manifest.get("source_sha256") != expected_source_hashes:
        raise InputError("export manifest has the wrong source hashes")
    output_hashes = manifest.get("output_sha256")
    if not isinstance(output_hashes, dict):
        raise InputError("export manifest has no output_sha256 mapping")
    for name in EXPECTED_OUTPUT_FILES:
        if name == manifest_path.name:
            continue
        expected = output_hashes.get(name)
        actual = sha256(export_dir / name)
        if expected != actual:
            raise InputError(
                f"export hash mismatch for {name}: {actual}; expected {expected}"
            )


def install_export(staging_dir: Path, output_dir: Path) -> None:
    """Atomically replace each validated data product in the destination."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in EXPECTED_OUTPUT_FILES:
        os.replace(staging_dir / name, output_dir / name)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("senscalc_root", type=Path, help="SensCalc Git checkout")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"export destination (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="verify the commit and source hashes without running Wolfram",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        verified = verify_sources(args.senscalc_root)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print(f"SensCalc {SENSCALC_TAG} production inputs verified")
    mx_system_ids = source_mx_system_ids(verified)
    for key, path in verified.items():
        suffix = f" (MX {mx_system_ids[key]})" if key in mx_system_ids else ""
        relative_path = path.relative_to(args.senscalc_root.resolve())
        print(f"  {key}: {relative_path}{suffix}")

    if mx_system_ids:
        systems = ", ".join(sorted(set(mx_system_ids.values())))
        print(
            "note: Wolfram MX is system-dependent; these inputs record "
            f"{systems}. The export manifest records the decoder system ID."
        )

    if args.check_only:
        return 0

    wolframscript = shutil.which("wolframscript")
    if wolframscript is None:
        print(
            "error: wolframscript is required to decode the upstream MX files",
            file=sys.stderr,
        )
        return 3

    output_dir = args.output_dir.resolve()
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    try:
        with tempfile.TemporaryDirectory(
            prefix=f".{output_dir.name}-", dir=output_dir.parent
        ) as staging_name:
            staging_dir = Path(staging_name)
            command = [
                wolframscript,
                "-file",
                str(EXPORTER),
                str(args.senscalc_root.resolve()),
                str(staging_dir),
            ]
            subprocess.run(command, check=True, env=wolfram_environment())
            validate_export(staging_dir)
            install_export(staging_dir, output_dir)
    except subprocess.CalledProcessError as exc:
        print(
            "error: Wolfram export failed. Activate Wolfram Engine locally "
            "with `wolframscript -activate` before retrying.",
            file=sys.stderr,
        )
        return exc.returncode or 4
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 5

    print(f"Exported SensCalc 2501 production inputs to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
