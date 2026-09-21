"""Durable writes and content hashes.

Every result a campaign keeps is written to a temporary file, flushed to
disk and renamed into place, so a crash never leaves a half-written CSV or
JSON where a complete one is expected.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path


def sha256_file(path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, payload, *, indent=2, sort_keys=True) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=indent, sort_keys=sort_keys, allow_nan=True)
        fh.write("\n")
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def atomic_csv(frame, path, **to_csv) -> None:
    """Write a pandas frame atomically (``index=False`` unless overridden)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    to_csv.setdefault("index", False)
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    frame.to_csv(tmp, **to_csv)
    with open(tmp, "rb") as fh:
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def tree_hash(entries) -> str:
    """Digest a file inventory (``path``, ``bytes``, ``sha256`` per entry)
    independently of where the tree lives on disk."""
    digest = hashlib.sha256()
    for entry in sorted(entries, key=lambda item: item["path"]):
        line = f"{entry['path']}\0{entry['bytes']}\0{entry['sha256']}\n"
        digest.update(line.encode())
    return digest.hexdigest()


def inventory(root, files) -> list[dict]:
    """``tree_hash`` entries for ``files`` relative to ``root``."""
    root = Path(root)
    out = []
    for f in files:
        f = Path(f)
        out.append({"path": str(f.relative_to(root)), "bytes": f.stat().st_size,
                    "sha256": sha256_file(f)})
    return out
