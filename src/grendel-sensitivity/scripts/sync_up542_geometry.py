"""Mirror the UP542 geometry from its single source into grendel/geometry/data/up542/, or check that the mirror is current.

    python scripts/sync_up542_geometry.py <clone of https://gitlab.cern.ch/grendel/grendel-geometry>           # sync
    python scripts/sync_up542_geometry.py --check <clone of https://gitlab.cern.ch/grendel/grendel-geometry>   # check only

The geometry is kept in one place, the GRENDEL geometry study (CERN GitLab, results/geometry/); this repository only
mirrors it. Sync checks the source against its manifest, replaces the mirror and records the source commit in
SOURCE.json. --check fetches the clone's remote (CERN login needed) and compares the latest results/geometry/ on its
main branch with the mirror, file by file (SHA-256); exit status 0 = up to date, 1 = the mirror is behind or differs.
CERN-internal: keep this repository private.
"""
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

DEST = Path(__file__).resolve().parents[1] / "grendel" / "geometry" / "data" / "up542"


def git(repo, *args):
    r = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"git {' '.join(args)} failed in {repo}: {r.stderr.strip()}")
    return r.stdout


def check(src_repo):
    git(src_repo, "fetch", "--quiet", "origin")
    ref = "origin/main"
    head = git(src_repo, "rev-parse", ref).strip()
    last_geo = git(src_repo, "log", "-1", "--format=%H %cs %s", ref, "--", "results/geometry").strip()
    remote = json.loads(git(src_repo, "show", f"{ref}:results/geometry/manifest.json"))
    local = json.loads((DEST / "manifest.json").read_text())
    source = json.loads((DEST / "SOURCE.json").read_text())
    intact = all(hashlib.sha256((DEST / f).read_bytes()).hexdigest() == h for f, h in local["files"].items())
    differ = sorted(f for f in set(remote["files"]) | set(local["files"])
                    if remote["files"].get(f) != local["files"].get(f))
    print(f"GitLab {ref} = {head[:9]}; last change to results/geometry: {last_geo[:9]} {last_geo[41:]}")
    print(f"mirror: version {local['version']}, copied from {source.get('commit', '?')[:9]} on {source.get('synced', '?')}; "
          f"intact: {intact}")
    if differ or not intact:
        print(f"NOT CURRENT: GitLab has version {remote['version']}; files that differ: {', '.join(differ) or '(none; mirror edited)'}")
        print("run: python scripts/sync_up542_geometry.py <clone>, then commit")
        return 1
    print(f"up to date: the mirror equals GitLab's results/geometry (version {remote['version']}, {len(local['files'])} files)")
    return 0


def sync(src_repo):
    src = src_repo / "results" / "geometry"
    man = json.loads((src / "manifest.json").read_text())
    for f, h in man["files"].items():
        if hashlib.sha256((src / f).read_bytes()).hexdigest() != h:
            sys.exit(f"source file {f} does not match its manifest: not syncing")
    commit = git(src_repo, "rev-parse", "HEAD").strip()
    DEST.mkdir(parents=True, exist_ok=True)
    for old in DEST.iterdir():
        if old.name not in man["files"] and old.name not in ("manifest.json", "SOURCE.json"):
            print("note: not in the source any more:", old.name)
    for f in list(man["files"]) + ["manifest.json"]:
        shutil.copy2(src / f, DEST / f)
    (DEST / "SOURCE.json").write_text(json.dumps(dict(
        mirror_of="https://gitlab.cern.ch/grendel/grendel-geometry", path="results/geometry/", commit=commit or "unknown",
        synced=str(date.today()), how="python scripts/sync_up542_geometry.py <path to a clone of grendel-geometry>; "
                                      "do not edit these files here"), indent=1) + "\n")
    print(f"mirrored UP542 geometry version {man['version']} ({len(man['files'])} files) from commit {commit[:9]}")
    return 0


def main():
    args = sys.argv[1:]
    if not args or args == ["--check"]:
        sys.exit(__doc__)
    if args[0] == "--check":
        sys.exit(check(Path(args[1]).resolve()))
    sys.exit(sync(Path(args[0]).resolve()))


if __name__ == "__main__":
    main()
