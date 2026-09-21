"""Tracked input tables, stored gzipped when large.

Loaders name the plain file; ``resolve_table`` returns whichever of
``<name>`` or ``<name>.gz`` exists and ``open_table`` opens it as text, so a
table can be compressed in the repository without touching its consumers.
"""
from __future__ import annotations

import gzip
import io
from pathlib import Path


def resolve_table(path) -> Path:
    path = Path(path)
    if path.exists():
        return path
    gz = path.with_name(path.name + ".gz")
    if gz.exists():
        return gz
    raise FileNotFoundError(f"neither {path} nor {gz.name} exists")


def open_table(path, mode: str = "rt", **kwargs):
    """Open ``path`` or its ``.gz`` twin. Text mode by default."""
    real = resolve_table(path)
    if real.suffix == ".gz":
        return gzip.open(real, mode, **kwargs)
    return open(real, mode, **kwargs)


def read_table_text(path) -> str:
    with open_table(path, "rt") as fh:
        return fh.read()


def read_table_bytes(path) -> bytes:
    with open_table(path, "rb") as fh:
        return fh.read()


def table_buffer(path) -> io.StringIO:
    """The table's text in a buffer, for readers that want a file object
    they can seek (``pandas.read_csv``, ``json.load``)."""
    return io.StringIO(read_table_text(path))
