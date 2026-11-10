"""Check the four-vector pools of step 02 byte for byte against the published runs.

The stored pool CSVs cannot be read here (macOS HFS+ compression: the copies on this disk are empty, see
campaigns/README.md "Retired and compressed"). Every scan records, per row, the SHA-256 over its input files
(grendel.scan.inputs_digest = sha256("mA_<m>.csv:<sha256>\\ntemplates_<m>.npz:<sha256>\\n")), and the template
bundles of those scans are stored. So a pool file is byte-identical to the published one iff the digest
rebuilt from our CSV and the stored template equals the recorded one. Two independent records:

  * v2.4 (the published table): bc10_v24_20261019/lxplus/out/bc10_<m>/run_metadata.json with the v2.4 templates
    of bc10_v24_20261019/templates_exhad (153 masses incl. the 8 tip masses);
  * v2.2: bc10_gkoz_cbs_20261005/analysis/run_metadata.json with the v2 templates of
    bc10/bc10_exhad_20260918/templates_exhad (the 143 grid masses).

    python validate_pools.py <work/alp/llp_4vectors> <report.txt>

The records and bundles are reference material ($GRENDEL_REF); without them the script prints SKIP and exits 0.
"""
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from ref import REF as REF_TREE, skip_unless  # noqa: E402  (pipelines/common/ref.py)

REF = REF_TREE / "grendel/campaigns"     # $GRENDEL_REF (optional)
V24 = REF / "bc10/bc10_v24_20261019"
GKOZ = REF / "bc10_gkoz_cbs_20261005"
V2_TEMPLATES = REF / "bc10/bc10_exhad_20260918/templates_exhad"


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def label(m):
    return f"{m:.3f}".replace(".", "p")


def digest(vec_sha, tpl_path, lab):
    h = hashlib.sha256()
    h.update(f"mA_{lab}.csv:{vec_sha}\n".encode())
    h.update(f"templates_{lab}.npz:{sha(tpl_path) if tpl_path.exists() else 'missing'}\n".encode())
    return h.hexdigest()


def main(vec_dir, report):
    skip_unless("BC10 pools vs the digests recorded by the published scans",
                V24 / "lxplus/points_all.txt", V24 / "lxplus/out", V24 / "templates_exhad",
                GKOZ / "analysis/run_metadata.json", V2_TEMPLATES)
    vec_dir = Path(vec_dir)
    gkoz = json.loads((GKOZ / "analysis/run_metadata.json").read_text())["inputs"]
    lines, n_ok, n_bad, n_missing = [], 0, 0, 0
    for text in (V24 / "lxplus/points_all.txt").read_text().split():
        m = float(text)
        lab, key = label(m), f"-/{m!r}"
        meta = V24 / f"lxplus/out/bc10_{text}/run_metadata.json"
        rec24 = json.loads(meta.read_text())["inputs"].get(key)
        vec = vec_dir / f"mA_{lab}.csv"
        if not vec.exists():                   # eta / eta' pole rows: no four-vectors, no templates
            tpl = V24 / f"templates_exhad/templates_{lab}.npz"
            if not tpl.exists() and digest("missing", tpl, lab) == rec24:
                lines.append(f"{m:<6g} pole row: no four-vectors here nor in the record (digest of 'missing' matches)")
            else:
                n_missing += 1
                lines.append(f"{m:<6g} MISSING")
            continue
        s = sha(vec)
        ok24 = digest(s, V24 / f"templates_exhad/templates_{lab}.npz", lab) == rec24
        ok22 = None
        if key in gkoz:
            ok22 = digest(s, V2_TEMPLATES / f"templates_{lab}.npz", lab) == gkoz[key]
        good = ok24 and ok22 is not False
        n_ok += good
        n_bad += not good
        lines.append(f"{m:<6g} sha256 {s[:16]}  v2.4 record {'match' if ok24 else 'DIFFER'}"
                     + ("" if ok22 is None else f"  v2.2 record {'match' if ok22 else 'DIFFER'}"))
    head = f"pools byte-identical to the published inputs: {n_ok} match, {n_bad} differ, {n_missing} missing"
    Path(report).write_text(head + "\n" + "\n".join(lines) + "\n")
    print(head)
    return 0 if n_bad == 0 and n_missing == 0 else 1


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
