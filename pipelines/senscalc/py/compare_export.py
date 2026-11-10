"""Compare a Wolfram export directory with a reference directory, file by file.

    python3 -I compare_export.py <ours dir> <reference dir> [--rtol 1e-12]

Every file of the reference directory (and every *.gz there, decompressed, against the uncompressed export) is
compared: first byte for byte; if different, numerically -- CSV cell by cell, JSON leaf by leaf (strings must be
equal, numbers within rtol). Files only in ours are listed. Exit 1 on any numerical difference.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import sys


def num_close(a, b, rtol):
    if a == b:
        return True
    if math.isnan(a) and math.isnan(b):
        return True
    return abs(a - b) <= rtol * max(abs(a), abs(b))


def cmp_json(a, b, rtol, path="$", out=None):
    out = [] if out is None else out
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            out.append(f"{path}: keys differ {sorted(set(a) ^ set(b))[:5]}")
        for k in sorted(set(a) & set(b)):  # sorted: the report lists only the first differences
            cmp_json(a[k], b[k], rtol, f"{path}.{k}", out)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(f"{path}: lengths {len(a)} vs {len(b)}")
        for i, (x, y) in enumerate(zip(a, b)):
            cmp_json(x, y, rtol, f"{path}[{i}]", out)
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        if not num_close(float(a), float(b), rtol):
            out.append(f"{path}: {a!r} vs {b!r}")
    elif a != b:
        out.append(f"{path}: {str(a)[:80]!r} vs {str(b)[:80]!r}")
    return out


def cmp_csv(ta, tb, rtol):
    ra, rb = list(csv.reader(io.StringIO(ta))), list(csv.reader(io.StringIO(tb)))
    if len(ra) != len(rb):
        return [f"rows {len(ra)} vs {len(rb)}"], 0.0
    bad, worst = [], 0.0
    for i, (x, y) in enumerate(zip(ra, rb)):
        if len(x) != len(y):
            bad.append(f"row {i}: {len(x)} vs {len(y)} cells")
            continue
        for j, (p, q) in enumerate(zip(x, y)):
            if p == q:
                continue
            try:
                fp, fq = float(p), float(q)
            except ValueError:
                bad.append(f"row {i} col {j}: {p!r} vs {q!r}")
                continue
            if fq != 0:
                worst = max(worst, abs(fp / fq - 1))
            if not num_close(fp, fq, rtol):
                bad.append(f"row {i} col {j}: {p} vs {q}")
    return bad, worst


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("ours")
    ap.add_argument("ref")
    ap.add_argument("--rtol", type=float, default=1e-12)
    a = ap.parse_args()
    status = 0
    for name in sorted(os.listdir(a.ref)):
        rp = os.path.join(a.ref, name)
        if not os.path.isfile(rp) or name in ("PROVENANCE.md", "LICENSE"):
            continue
        raw_ref = open(rp, "rb").read()
        op = os.path.join(a.ours, name)
        line = f"{name:48s}"
        if os.path.isfile(op):
            same = hashlib.sha256(open(op, "rb").read()).hexdigest() == hashlib.sha256(raw_ref).hexdigest()
            line += f" file {'IDENTICAL' if same else 'differs'};"
        if name.endswith(".gz"):
            ref_txt = gzip.decompress(raw_ref)
            base = name[:-3]
            op = os.path.join(a.ours, base)
        else:
            ref_txt, base = raw_ref, name
        if not os.path.isfile(op):
            print(f"{line} MISSING in ours ({base})")
            status = 1
            continue
        ours_txt = open(op, "rb").read()
        if ours_txt == ref_txt:
            print(f"{line} content IDENTICAL (sha256 {hashlib.sha256(ref_txt).hexdigest()[:16]})")
            continue
        if base.endswith(".csv"):
            bad, worst = cmp_csv(ours_txt.decode(), ref_txt.decode(), a.rtol)
            print(f"{line} content differs: CSV max |ours/ref-1| = {worst:.2e}; {len(bad)} cells beyond rtol {a.rtol:g}")
        else:
            bad = cmp_json(json.loads(ours_txt), json.loads(ref_txt), a.rtol)
            print(f"{line} content differs: JSON {len(bad)} leaves beyond rtol {a.rtol:g} (strings exact)")
        for b in bad[:10]:
            print("    ", b)
        if bad:
            status = 1
    extra = sorted(set(os.listdir(a.ours)) - set(os.listdir(a.ref)) - {n[:-3] for n in os.listdir(a.ref) if n.endswith(".gz")})
    if extra:
        print("only in ours:", ", ".join(extra))
    print("RESULT:", "PASS" if status == 0 else "FAIL")
    return status


if __name__ == "__main__":
    sys.exit(main())
