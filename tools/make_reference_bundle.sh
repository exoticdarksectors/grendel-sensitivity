#!/usr/bin/env bash
# tools/make_reference_bundle.sh DEST -- pack the reference files the validations read (reference/MANIFEST.tsv) from
# $GRENDEL_REF into one tar archive per pipeline, DEST/grendel-sensitivity-reference-<pipeline>.tar, for a data record
# (e.g. Zenodo). Unpacked into one folder, the archives rebuild the layout of the reference tree; point GRENDEL_REF
# at that folder (reference/README.md). The files of reference/not_redistributable.txt (given to the authors
# privately) are left out.
#
#     source env/activate.sh && bash tools/make_reference_bundle.sh /path/with/space [pipeline ...]
#
# Without pipeline names, all of them (36 GB in all, reference/README.md).
#
# Nothing is uploaded. The archives are checked against the manifest after writing (size and SHA-256 per file).
set -euo pipefail
REPRO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${1:?usage: make_reference_bundle.sh DEST [pipeline ...]}"
shift
: "${GRENDEL_REF:?source env/activate.sh first}"
[ "$GRENDEL_REF" != none ] && [ -d "$GRENDEL_REF" ] || { echo "GRENDEL_REF is not a reference tree" >&2; exit 2; }
M="$REPRO/reference/MANIFEST.tsv"
X="$REPRO/reference/not_redistributable.txt"
mkdir -p "$DEST"
PIPES=("$@")
[ ${#PIPES[@]} -gt 0 ] || mapfile -t PIPES < <(tail -n +2 "$M" | cut -f1 | sort -u)
for pipe in "${PIPES[@]}"; do
  out="$DEST/grendel-sensitivity-reference-$pipe.tar"
  tail -n +2 "$M" | awk -F'\t' -v p="$pipe" '$1 == p { print $2 }' | grep -vxF -f <(grep -v '^#' "$X") \
    > "$DEST/.$pipe.list" || true
  [ -s "$DEST/.$pipe.list" ] || { rm -f "$DEST/.$pipe.list"; continue; }
  tar -C "$GRENDEL_REF" --owner=0 --group=0 --numeric-owner -cf "$out" -T "$DEST/.$pipe.list"
  # check: every member has the size and SHA-256 of the manifest
  python3 - "$out" "$M" "$pipe" "$X" <<'EOF'
import hashlib, sys, tarfile
out, manifest, pipe, excluded = sys.argv[1:]
skip = {l.strip() for l in open(excluded) if l.strip() and not l.startswith("#")}
want = {}
for line in open(manifest).read().splitlines()[1:]:
    p, rel, size, digest = line.split("\t")
    if p == pipe and rel not in skip:
        want[rel] = (int(size), digest)
seen = 0
with tarfile.open(out) as t:
    for m in t:
        if not m.isfile():
            continue
        h = hashlib.sha256()
        f = t.extractfile(m)
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
        assert want.get(m.name) == (m.size, h.hexdigest()), m.name
        seen += 1
assert seen == len(want), (seen, len(want))
print(f"{out}: {seen} files, checked")
EOF
  rm -f "$DEST/.$pipe.list"
done
(cd "$DEST" && sha256sum grendel-sensitivity-reference-*.tar > SHA256SUMS && du -ch grendel-sensitivity-reference-*.tar | tail -1)
