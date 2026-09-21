"""Rebuild ``variation_manifest.json`` from the campaign grids on disk.

    python -m grendel.production.fonll_grids.snapshot [--out-dir DIR]

``generate --campaign`` writes its manifest only when the whole run
finishes. This scans the per-variation ``.dat`` grids already written and
emits an equivalent manifest for whatever is complete now, so a coherent
set can be combined and banked at any time (hourly, or after an
interrupted run). Every field comes from the grids' own headers; only the
coherent ``as_01180`` campaign grids are included (alpha_s companions,
envelopes and diagnostic outputs are skipped).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from . import generate as gen

PDF_SET = "NNPDF40_nlo_as_01180"
PDF_BASE_ID = 331700
STABLE_AGE_S = 5.0  # skip files still being written


def parse_header(path: Path) -> dict | None:
    """The variation fields of a grid header; None if not a campaign grid."""
    hdr: dict[str, str] = {}
    rows = 0
    with path.open() as f:
        for line in f:
            if line.startswith("#"):
                body = line[1:].strip()
                if ":" in body:
                    k, v = body.split(":", 1)
                    hdr[k.strip()] = v.strip()
            elif line.strip():
                rows += 1
    if hdr.get("pdf") != PDF_SET or "variation_kind" not in hdr:
        return None
    # "scale: mu0=sqrt(m^2+pT^2); ffact(muF)=1, fren(muR)=1"
    muF = muR = 1.0
    for tok in hdr.get("scale", "").replace(";", ",").split(","):
        if "ffact(muF)=" in tok:
            muF = float(tok.split("=", 1)[1])
        elif "fren(muR)=" in tok:
            muR = float(tok.split("=", 1)[1])
    return {"quark": hdr["quark"], "variation_kind": hdr["variation_kind"],
            "variation_tag": hdr.get("variation_tag", ""), "muR": muR, "muF": muF,
            "lhapdf_id": int(hdr["lhapdf_id"]), "lhapdf_member": int(hdr["lhapdf_member"]),
            "heavy_quark_mass_GeV": float(hdr["heavy_quark_mass_GeV"]), "path": str(path),
            "sha256": gen.sha256_file(path), "rows": rows,
            "trapezoid_integral_pb": float(hdr["trapezoid_integral_pb_y-3to3_pt0to50"])}


def snapshot(ws: gen.Workspace) -> Path:
    now = time.time()
    entries = []
    for path in sorted(ws.out.glob("*.dat")):
        name = path.name
        # Campaign grids end exactly "_bottom.dat"/"_charm.dat" with no extra
        # dotted suffix (diagnostic variants have ".something.dat").
        if name.count(".") != 1 or not (name.endswith("_bottom.dat") or name.endswith("_charm.dat")):
            continue
        if now - path.stat().st_mtime < STABLE_AGE_S:
            continue
        entry = parse_header(path)
        if entry is not None:
            entries.append(entry)
    entries.sort(key=lambda e: (e["quark"], e["variation_tag"]))
    members = sorted({e["lhapdf_member"] for e in entries})
    feeddown_path = ws.out / "charm_feeddown_calibration.json"
    feeddown = json.loads(feeddown_path.read_text()) if feeddown_path.exists() else None
    manifest = gen.manifest_document(ws, entries, PDF_SET, PDF_BASE_ID, members, gen.SCALE_POINTS_7,
                                     gen.MASS_VARIATIONS, feeddown,
                                     note="on-disk snapshot (partial campaign allowed)")
    out = ws.out / "variation_manifest.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n")
    by_q: dict[str, list[int]] = {}
    for e in entries:
        by_q.setdefault(e["quark"], []).append(e["lhapdf_member"])
    print(f"snapshot: {len(entries)} campaign grids -> {out} | "
          + "; ".join(f"{q}: {len(m)} grids, replicas<= {max([x for x in m if x >= 1], default=0)}"
                      for q, m in by_q.items()))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    gen.add_workspace_options(ap)
    args = ap.parse_args(argv)
    snapshot(gen.workspace_from(args))
    return 0


if __name__ == "__main__":
    sys.exit(main())
