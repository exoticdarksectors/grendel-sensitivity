"""Helpers of 03_templates.sh (run with exHad's venv python; standard library only).

    merge_manifests.py check <templates_exhad> <templates_jobs>
        remove template bundles that do not open (a job killed while writing), so --skip-existing never keeps one;
    merge_manifests.py merge <templates_exhad> <templates_jobs>
        move the bundles of the job folders into templates_exhad and merge the jobs' MANIFEST.json records (which
        carry the retry logs) and generation runs into templates_exhad/MANIFEST.json; the driver's final pass then
        writes the complete manifest from it.
"""
import json
import sys
import zipfile
from pathlib import Path


def check(out, jobs):
    for p in sorted(out.glob("templates_*.npz")) + sorted(jobs.glob("job*/templates_*.npz")):
        try:
            with zipfile.ZipFile(p) as z:
                bad = z.testzip()
        except Exception as err:   # noqa: BLE001 -- any failure to open means a truncated file
            bad = err
        if bad is not None:
            print(f"removing unreadable {p}: {bad}")
            p.unlink()


def merge(out, jobs):
    man = out / "MANIFEST.json"
    base = json.loads(man.read_text()) if man.exists() else {}
    points = {p["file"]: p for p in base.get("points", []) if "file" in p}
    runs = list(base.get("generation_runs", []))
    job_manifests = sorted(jobs.glob("job*/MANIFEST.json"))
    for jm in job_manifests:
        m = json.loads(jm.read_text())
        if not base:
            base = m
        for p in m.get("points", []):
            if "file" in p and (jm.parent / p["file"]).exists():
                points[p["file"]] = p
        runs += [dict(r, job=jm.parent.name) for r in m.get("generation_runs", [])[-1:]]
    moved = 0
    for b in sorted(jobs.glob("job*/templates_*.npz")):
        b.replace(out / b.name)
        moved += 1
    for jm in job_manifests:
        jm.replace(jm.with_name("MANIFEST.merged.json"))
    if base:
        base["points"] = sorted(points.values(), key=lambda p: p["mass_GeV"])
        base["generation_runs"] = runs
        man.write_text(json.dumps(base, indent=2))
    print(f"moved {moved} bundles; manifest records {len(points)}")


if __name__ == "__main__":
    action, out, jobs = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    {"check": check, "merge": merge}[action](out, jobs)
