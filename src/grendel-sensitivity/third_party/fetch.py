#!/usr/bin/env python3
"""Fetch the external software the pipeline drives, at the pinned revisions.

    python third_party/fetch.py --list
    python third_party/fetch.py exHad HNLCalc SM_HeavyN_CKM_AllMasses_LO
    python third_party/fetch.py --all
    python third_party/fetch.py --verify

Nothing here is redistributed: each package is checked out from its own
upstream into ``third_party/<name>`` (or ``GRENDEL_THIRD_PARTY_DIR``) at the
commit or archive hash recorded in ``PINS.json``, which is what the
published results were produced with. Git sources are fetched shallowly at
the pinned commit (``sparse`` entries fetch only the listed paths); archive
sources are downloaded, hashed and unpacked. Every git checkout gets an
``UPSTREAM_COMMIT.txt`` with the commit, so a build can record it.

Fetching is not building. ``--list`` prints, per package, the step that
follows (FONLL patch + build, exHad's ``tools/setup.py``, Pythia's
``configure && make`` ...); THIRD_PARTY.md at the repository root carries
the licence notices and what each package is used for.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
PINS_PATH = HERE / "PINS.json"


def load_pins() -> dict:
    return json.loads(PINS_PATH.read_text())["packages"]


def dest_root() -> Path:
    return Path(os.environ.get("GRENDEL_THIRD_PARTY_DIR", HERE)).expanduser()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd: list[str], cwd: Path | None = None) -> str:
    return subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True).stdout.strip()


# ------------------------------------------------------------------ fetch --

def fetch_git(name: str, pin: dict, dest: Path) -> None:
    dest.mkdir(parents=True)
    run(["git", "init", "-q"], dest)
    run(["git", "remote", "add", "origin", pin["url"]], dest)
    fetch = ["git", "fetch", "-q", "--depth", "1"]
    sparse = pin.get("sparse")
    if sparse:
        fetch.append("--filter=blob:none")
    run([*fetch, "origin", pin["commit"]], dest)
    if sparse:
        run(["git", "sparse-checkout", "init", "--cone"], dest)
        run(["git", "sparse-checkout", "set", *sparse], dest)
    run(["git", "checkout", "-q", "FETCH_HEAD"], dest)
    head = run(["git", "rev-parse", "HEAD"], dest)
    if head != pin["commit"]:
        raise RuntimeError(f"{name}: checked out {head}, pinned {pin['commit']}")
    (dest / "UPSTREAM_COMMIT.txt").write_text(head + "\n")


def fetch_archive(name: str, pin: dict, dest: Path) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        archive = Path(tmp) / "archive"
        with urllib.request.urlopen(pin["url"], timeout=120) as response, archive.open("wb") as out:
            shutil.copyfileobj(response, out)
        digest = sha256_file(archive)
        if digest != pin["sha256"]:
            raise RuntimeError(f"{name}: archive sha256 {digest} != pinned {pin['sha256']}")
        with tarfile.open(archive) as tar:
            tar.extractall(tmp, filter="data")
        top = Path(tmp) / pin["archive_dir"]
        if not top.is_dir():
            raise RuntimeError(f"{name}: {pin['archive_dir']} not found in the archive")
        shutil.move(str(top), str(dest))
    (dest / "UPSTREAM_ARCHIVE.txt").write_text(f"{pin['url']}\nsha256 {pin['sha256']}\n")


def fetch(name: str, pin: dict, root: Path, force: bool) -> Path:
    dest = root / pin.get("dest", name)
    if dest.exists():
        if not force:
            print(f"{name}: already present at {dest} (use --force to refetch)")
            return dest
        shutil.rmtree(dest)
    print(f"{name}: fetching {pin['url']} @ {pin.get('commit') or pin.get('sha256')[:12]} -> {dest}")
    if pin["kind"] == "git":
        fetch_git(name, pin, dest)
    elif pin["kind"] == "archive":
        fetch_archive(name, pin, dest)
    else:
        raise ValueError(f"{name}: unknown kind {pin['kind']!r}")
    if pin.get("next"):
        print(f"  next: {pin['next']}")
    return dest


# ----------------------------------------------------------------- verify --

def verify(name: str, pin: dict, root: Path) -> tuple[str, str]:
    """(state, detail): ``ok``, ``missing`` or ``mismatch``."""
    dest = root / pin.get("dest", name)
    if not dest.is_dir():
        return "missing", str(dest)
    if pin["kind"] == "git":
        marker = dest / "UPSTREAM_COMMIT.txt"
        head = marker.read_text().strip() if marker.is_file() else None
        if head is None and (dest / ".git").exists():
            try:
                head = run(["git", "rev-parse", "HEAD"], dest)
            except subprocess.CalledProcessError:
                head = None
        if head is None:
            return "mismatch", f"{dest}: no UPSTREAM_COMMIT.txt and no git history"
        return ("ok" if head == pin["commit"] else "mismatch"), f"{dest} @ {head[:12]}"
    marker = dest / "UPSTREAM_ARCHIVE.txt"
    if marker.is_file() and pin["sha256"] in marker.read_text():
        return "ok", f"{dest} (archive {pin['sha256'][:12]})"
    return "mismatch", f"{dest}: no UPSTREAM_ARCHIVE.txt recording {pin['sha256'][:12]}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("names", nargs="*", help="packages to fetch (see --list)")
    ap.add_argument("--all", action="store_true", help="fetch every package")
    ap.add_argument("--list", action="store_true", help="list the packages, pins and next steps")
    ap.add_argument("--verify", action="store_true", help="report which packages are present at their pins")
    ap.add_argument("--force", action="store_true", help="refetch a package that is already present")
    ap.add_argument("--dest", type=Path, default=None, help=f"checkout root (default: {dest_root()})")
    args = ap.parse_args(argv)
    pins = load_pins()
    root = (args.dest or dest_root()).expanduser()

    if args.list:
        for name, pin in pins.items():
            pinned = pin.get("commit") or pin.get("sha256")
            print(f"{name}\n  {pin['what']}\n  {pin['url']}\n  {pin['kind']} {pinned}  licence: {pin['license']}")
            if pin.get("next"):
                print(f"  next: {pin['next']}")
        return 0
    if args.verify:
        worst = 0
        for name, pin in pins.items():
            state, detail = verify(name, pin, root)
            print(f"{state:9s} {name:32s} {detail}")
            worst = max(worst, {"ok": 0, "missing": 1, "mismatch": 2}[state])
        return worst
    names = list(pins) if args.all else args.names
    if not names:
        ap.error("name at least one package, or pass --all / --list / --verify")
    unknown = [n for n in names if n not in pins]
    if unknown:
        ap.error(f"unknown package(s): {', '.join(unknown)}; see --list")
    root.mkdir(parents=True, exist_ok=True)
    for name in names:
        fetch(name, pins[name], root, args.force)
    return 0


if __name__ == "__main__":
    sys.exit(main())
