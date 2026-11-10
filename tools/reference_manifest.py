"""Write reference/MANIFEST.tsv: size and SHA-256 of every file of the reference tree ($GRENDEL_REF) that a
validation step reads. The lists of those files, reference/files/<pipeline>.txt (paths relative to $GRENDEL_REF),
were recorded by running the validations under strace.

    source env/activate.sh && python tools/reference_manifest.py          # needs the reference tree
    python tools/reference_manifest.py --check [pipeline ...]             # verify a tree against the manifest

Columns: pipeline, path, size in bytes, sha256. A path read by several pipelines is listed once per pipeline.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
from pathlib import Path

REPRO = Path(__file__).resolve().parents[1]
LISTS = REPRO / "reference/files"
MANIFEST = REPRO / "reference/MANIFEST.tsv"
PRIVATE = REPRO / "reference/not_redistributable.txt"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def entries():
    for lst in sorted(LISTS.glob("*.txt")):
        for line in lst.read_text().splitlines():
            if line.strip() and not line.startswith("#"):
                yield lst.stem, line.strip()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="verify $GRENDEL_REF against the manifest")
    ap.add_argument("pipelines", nargs="*", help="with --check: only these pipelines (default: all)")
    a = ap.parse_args()
    root = os.environ.get("GRENDEL_REF", "")
    if root in ("", "none") or not Path(root).is_dir():
        print("GRENDEL_REF is not set to a reference tree", file=sys.stderr)
        return 2
    root = Path(root)
    if a.check:
        private = {l.strip() for l in PRIVATE.read_text().splitlines() if l.strip() and not l.startswith("#")}
        bad = n = absent = 0
        for line in MANIFEST.read_text().splitlines()[1:]:
            pipe, rel, size, digest = line.split("\t")
            if a.pipelines and pipe not in a.pipelines:
                continue
            n += 1
            p = root / rel
            if not p.is_file() and rel in private:
                print(f"absent   {rel} (not redistributable; its comparisons are skipped)")
                absent += 1
            elif not p.is_file():
                print(f"missing  {rel}")
                bad += 1
            elif p.stat().st_size != int(size) or sha256(p) != digest:
                print(f"differs  {rel}")
                bad += 1
        print(f"{n - bad - absent} of {n} files as in the manifest" + (f", {absent} not redistributable" if absent else ""))
        return 1 if bad else 0
    rows, cache = [], {}
    for pipe, rel in entries():
        p = root / rel
        if not p.is_file():
            print(f"not a file in {root}: {rel}", file=sys.stderr)
            return 1
        if rel not in cache:
            cache[rel] = (p.stat().st_size, sha256(p))
        rows.append((pipe, rel, *cache[rel]))
    with MANIFEST.open("w") as f:
        f.write("pipeline\tpath\tsize\tsha256\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    total = sum(s for s, _ in cache.values())
    empty = sum(1 for s, _ in cache.values() if s == 0)
    print(f"{MANIFEST.relative_to(REPRO)}: {len(rows)} rows, {len(cache)} distinct files, {total / 1e9:.2f} GB"
          f" ({empty} empty)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
