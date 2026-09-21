#!/usr/bin/env python3
"""Refuse absolute paths and machine-specific references in the tree.

    python scripts/check_no_private_refs.py            # the working tree
    python scripts/check_no_private_refs.py --history  # every commit message too
    python scripts/check_no_private_refs.py PATH ...   # specific files

Every input the pipeline needs is resolved through ``grendel.io.paths``
(environment variables with documented defaults) or fetched by
``third_party/fetch.py``; a hard-coded ``/Volumes/...``, a conda prefix or a
reference to a private checkout would only work on one machine. This scan
fails CI when one slips in. Compressed data files are scanned decompressed.
"""
from __future__ import annotations

import argparse
import gzip
import re
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

PATTERNS = [
    (r"/Volumes/", "absolute macOS volume path"),
    (r"/Users/[a-z]", "absolute home path"),
    (r"/home/[a-z]", "absolute home path"),
    (r"conda/envs/", "conda environment path"),
    (r"[Cc]o-[Aa]uthored-[Bb]y:", "commit trailer"),
]

# Files whose content legitimately carries such strings.
ALLOW = {
    "scripts/check_no_private_refs.py",
}
# Binary or generated files that are not text.
SKIP_SUFFIXES = {".png", ".pdf", ".npz", ".pyc", ".so", ".dylib", ".ico", ".jpg", ".pkl"}


def tracked_files() -> list[Path]:
    out = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, check=True, capture_output=True).stdout
    return [ROOT / p for p in out.decode().split("\0") if p]


def contents(path: Path) -> list[tuple[str, bytes]]:
    """(label, bytes) for a file, decompressing .gz and .xlsx/.zip members."""
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as fh:
            return [(f"{path}", fh.read())]
    if path.suffix in (".xlsx", ".zip"):
        with zipfile.ZipFile(path) as zf:
            return [(f"{path}!{name}", zf.read(name)) for name in zf.namelist()]
    return [(str(path), path.read_bytes())]


def scan_blob(label: str, data: bytes) -> list[str]:
    hits = []
    text = data.decode("utf-8", errors="replace")
    for pattern, what in PATTERNS:
        for match in re.finditer(pattern, text):
            line = text.count("\n", 0, match.start()) + 1
            snippet = text[max(0, match.start() - 30): match.end() + 40].replace("\n", " ")
            hits.append(f"{label}:{line}: {what}: ...{snippet}...")
    return hits


def scan_files(paths: list[Path]) -> list[str]:
    hits = []
    for path in paths:
        rel = path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else str(path)
        if rel in ALLOW or path.suffix in SKIP_SUFFIXES or not path.is_file():
            continue
        for label, data in contents(path):
            hits.extend(scan_blob(label.replace(str(ROOT) + "/", ""), data))
    return hits


def scan_history() -> list[str]:
    marker = "----commit----"
    log = subprocess.run(["git", "log", "--all", f"--format={marker}%n%H%n%B"], cwd=ROOT, check=True,
                         capture_output=True).stdout.decode()
    hits = []
    for entry in log.split(marker):
        entry = entry.strip()
        if not entry:
            continue
        sha, _, body = entry.partition("\n")
        hits.extend(f"commit {sha[:12]}: {h}" for h in scan_blob("message", body.encode()))
    return hits


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="*", type=Path, help="files to scan (default: every tracked file)")
    ap.add_argument("--history", action="store_true", help="also scan every commit message")
    args = ap.parse_args(argv)
    paths = [p.resolve() for p in args.paths] if args.paths else tracked_files()
    hits = scan_files(paths)
    if args.history:
        hits.extend(scan_history())
    for hit in hits:
        print(hit)
    print(f"{len(paths)} files scanned, {len(hits)} problems")
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main())
