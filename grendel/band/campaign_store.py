"""Bookkeeping shared by the uncertainty campaigns.

A campaign runs hundreds of independent variations for days, so every
retained product is written atomically, every run directory is locked
against a second worker, and completion is recorded in markers that pin
the configuration, the code and the content hashes of what was produced.
"""
from __future__ import annotations

import fcntl
import hashlib
import json
import os
import platform
import shutil
import subprocess
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from ..io.atomic import sha256_file, tree_hash
from ..io.paths import repo_root


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def jsonable(value):
    """Plain JSON types for numpy scalars, paths and nested containers."""
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def config_sha256(payload: dict) -> str:
    stable = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(stable).hexdigest()


def stable_seed(base_seed: int, *parts) -> int:
    """A 32-bit seed derived from the base seed and any labels."""
    message = ":".join([str(base_seed), *map(str, parts)]).encode()
    return int.from_bytes(hashlib.sha256(message).digest()[:8], "big") % (2**32)


def git_head(cwd=None) -> str:
    proc = subprocess.run(["git", "rev-parse", "HEAD"], cwd=cwd or repo_root(),
                          text=True, capture_output=True, check=False)
    return proc.stdout.strip() if proc.returncode == 0 else "unknown"


def git_state(cwd=None) -> dict:
    """Commit plus a digest of the uncommitted tracked diff."""
    cwd = cwd or repo_root()
    commit = git_head(cwd)
    proc = subprocess.run(["git", "diff", "--binary", "HEAD"], cwd=cwd,
                          capture_output=True, check=False)
    diff = proc.stdout if proc.returncode == 0 else b""
    return {"commit": commit,
            "tracked_diff_sha256": hashlib.sha256(diff).hexdigest(),
            "tracked_tree_clean": not bool(diff)}


def code_hashes(files) -> dict[str, str]:
    """sha256 of repository files (paths relative to the repository root)."""
    root = repo_root()
    return {rel: sha256_file(root / rel) for rel in files}


def validate_recorded_code_state(commit: str, hashes: dict, expected_files) -> str:
    """Prove that retained producer hashes match the files in the named commit."""
    root = repo_root()
    resolved = subprocess.check_output(["git", "rev-parse", f"{commit}^{{commit}}"],
                                       cwd=root, text=True).strip()
    if set(hashes) != set(expected_files):
        raise RuntimeError("recorded producer code inventory is incomplete")
    for relative in expected_files:
        recorded = str(hashes[relative])
        if len(recorded) != 64 or any(c not in "0123456789abcdef" for c in recorded):
            raise RuntimeError(f"invalid recorded code checksum for {relative}")
        content = subprocess.check_output(["git", "show", f"{resolved}:{relative}"], cwd=root)
        if hashlib.sha256(content).hexdigest() != recorded:
            raise RuntimeError(f"recorded producer checksum does not match {resolved}:{relative}")
    return resolved


@contextmanager
def variation_lock(run_dir: Path):
    """An exclusive, non-blocking lock on a variation's run directory."""
    run_dir.mkdir(parents=True, exist_ok=True)
    lock_path = run_dir / ".lock"
    with open(lock_path, "a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"variation already running: {run_dir.name}") from exc
        lock.seek(0)
        lock.truncate()
        lock.write(f"pid={os.getpid()} host={platform.node()}\n")
        lock.flush()
        try:
            yield
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def append_log(run_dir: Path, message: str) -> None:
    """One durable timestamped line in the variation's retained log."""
    run_dir.mkdir(parents=True, exist_ok=True)
    with open(run_dir / "run.log", "a") as fh:
        fh.write(f"{utc_now()} {message}\n")
        fh.flush()
        os.fsync(fh.fileno())


def file_record(path: Path, run_dir: Path, known_sha=None) -> dict:
    return {"path": str(path.relative_to(run_dir)), "bytes": path.stat().st_size,
            "sha256": known_sha or sha256_file(path)}


def stage_record(entries) -> dict:
    entries = sorted(entries, key=lambda item: item["path"])
    return {"n_files": len(entries), "bytes": int(sum(e["bytes"] for e in entries)),
            "tree_sha256": tree_hash(entries), "files": entries}


def read_json(path: Path):
    """The JSON at ``path`` or None when absent or unreadable."""
    if not Path(path).exists():
        return None
    try:
        return json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError):
        return None


def remove_generated_tree(path: Path, run_dir: Path) -> None:
    """Delete a reclaimable tree, refusing anything outside the run directory."""
    path = Path(path).resolve()
    run_dir = Path(run_dir).resolve()
    if path.parent != run_dir and path.parent.parent != run_dir:
        raise RuntimeError(f"refusing to compact path outside run directory: {path}")
    if path.exists():
        shutil.rmtree(path)


def tree_usage(path: Path) -> dict:
    files = [item for item in Path(path).rglob("*") if item.is_file()]
    return {"path": str(path), "n_files": len(files),
            "bytes": sum(item.stat().st_size for item in files)}


def software() -> dict:
    import pandas
    import sys
    return {"python": sys.version.split()[0], "numpy": np.__version__,
            "pandas": pandas.__version__, "platform": platform.platform()}
