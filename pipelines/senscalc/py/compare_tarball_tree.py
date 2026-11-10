"""Check that a GitHub tarball holds exactly the tree of a git commit (git blob hashes, file by file).

    python3 -I compare_tarball_tree.py <tarball.tar.gz> <git ls-tree -r output>

Prints the number of files compared and every mismatch; exit status 1 on any difference.
"""
import hashlib
import sys
import tarfile


def main() -> int:
    tgz, lstree = sys.argv[1], sys.argv[2]
    tree = {}
    for line in open(lstree, encoding="utf-8"):
        meta, path = line.rstrip("\n").split("\t", 1)
        mode, kind, sha = meta.split()
        if path.startswith('"'):  # git quotes unusual paths
            path = path[1:-1].encode("latin-1").decode("unicode_escape").encode("latin-1").decode("utf-8")
        if kind == "blob":
            tree[path] = sha
    seen, bad = 0, []
    with tarfile.open(tgz, "r:gz") as tf:
        for m in tf:
            if not (m.isfile() or m.issym()):
                continue
            rel = m.name.split("/", 1)[1]
            data = m.linkname.encode() if m.issym() else tf.extractfile(m).read()
            h = hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()
            seen += 1
            if tree.get(rel) != h:
                bad.append((rel, tree.get(rel), h))
            tree.pop(rel, None)
    print(f"tarball files compared: {seen}; mismatches: {len(bad)}; tree files missing from tarball: {len(tree)}")
    for b in bad[:20]:
        print("  MISMATCH", b)
    for p in list(tree)[:20]:
        print("  ONLY IN TREE", p)
    return 1 if bad or tree else 0


if __name__ == "__main__":
    sys.exit(main())
