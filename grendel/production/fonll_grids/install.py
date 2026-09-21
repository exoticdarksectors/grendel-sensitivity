"""Apply the patches to a FONLL v1.3.3 tree and build the two executables.

    python -m grendel.production.fonll_grids.install --fonll DIR [--build]

``DIR`` is an unpacked FONLL v1.3.3 source tree (``main/``, ``misc1/``,
``Linux/`` ...), either the authors' tarball or the ``alisw/fonll`` mirror
at the commit pinned in ``patches/PINS.json``. Each file's SHA-256 is
checked before patching (it must be pristine) and after (it must match
the patched hash the grids were produced with); a tree that is already
patched is recognised and left alone.

``--build`` then runs the FONLL Makefile for ``fonllgridlha`` (the grid
driver linked against LHAPDF, needs ``lhapdf-config`` on PATH) and
``fragmfonll`` (the fragmentation convolution) in ``DIR/Linux``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

from ...io.paths import fonll_dir

HERE = Path(__file__).resolve().parent
PATCHES = HERE / "patches"
PINS = json.loads((PATCHES / "PINS.json").read_text())
EXECUTABLES = ("fonllgridlha", "fragmfonll")
MAKE = ["make", "-f", "../misc1/Makefile",
        "VPATH=../misc1:../main:../hdmassive:../hdresummed:../phmassive:../phresummed:../common"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def status(tree: Path) -> dict[str, str]:
    """``pristine`` / ``patched`` / ``unknown`` / ``missing`` per patched file."""
    out = {}
    for rel, pin in PINS["files"].items():
        path = tree / rel
        if not path.is_file():
            out[rel] = "missing"
            continue
        digest = sha256(path)
        out[rel] = ("pristine" if digest == pin["sha256_pristine"]
                    else "patched" if digest == pin["sha256_patched"] else "unknown")
    return out


def apply_patches(tree: Path) -> None:
    state = status(tree)
    for rel, s in state.items():
        if s in ("missing", "unknown"):
            raise SystemExit(f"{tree / rel}: {s}; expected a pristine FONLL {PINS['fonll']['version']} file "
                             f"(mirror commit {PINS['fonll']['mirror_commit'][:12]})")
    if all(s == "patched" for s in state.values()):
        print("already patched")
        return
    for rel, pin in PINS["files"].items():
        if state[rel] == "patched":
            continue
        diff = PATCHES / pin["diff"]
        subprocess.run(["patch", "-p1", "-s", "-d", str(tree)], input=diff.read_bytes(), check=True)
        digest = sha256(tree / rel)
        if digest != pin["sha256_patched"]:
            raise SystemExit(f"{rel}: patched file hash {digest} != pinned {pin['sha256_patched']}")
        print(f"patched {rel}")


def build(tree: Path) -> None:
    linux = tree / "Linux"
    linux.mkdir(exist_ok=True)
    subprocess.run([*MAKE, *EXECUTABLES], cwd=linux, check=True)
    for exe in EXECUTABLES:
        if not (linux / exe).is_file():
            raise SystemExit(f"build did not produce {linux / exe}")
    print(f"built {', '.join(EXECUTABLES)} in {linux}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fonll", type=Path, default=None, help=f"FONLL source tree (default: {fonll_dir()})")
    ap.add_argument("--build", action="store_true", help="also build fonllgridlha and fragmfonll")
    ap.add_argument("--status", action="store_true", help="report the patch state of each file and exit")
    args = ap.parse_args(argv)
    tree = (args.fonll or fonll_dir()).expanduser().resolve()
    if not tree.is_dir():
        raise SystemExit(f"FONLL tree not found: {tree}")
    if args.status:
        for rel, s in status(tree).items():
            print(f"{s:9s} {rel}")
        return 0
    apply_patches(tree)
    if args.build:
        build(tree)
    return 0


if __name__ == "__main__":
    sys.exit(main())
