"""Step 05: compare the regenerated FONLL outputs with the tracked grids and the stored original (macOS) run.

    source env/activate.sh && source pipelines/fonll/env.sh && python pipelines/fonll/05_compare.py

Checks whatever exists (it can run after step 02 alone), writes work/fonll/validation/compare.json and prints a
summary (also work/fonll/validation/summary.txt):

1. meson grids (step 04) vs the tracked src/.../data/fonll/central grids: full-file SHA-256, header differences,
   SHA-256 of the data rows alone, identical rows, per-bin relative difference (max / 99th percentile / median,
   the pT = 0 rows must be exactly 0 in both), and the integrated cross section (trapezoid over the grid, as
   grendel.production.fonll.fonll_parser.get_sigma_total), in total and in pT / |y| slices;
2. raw quark grids (fonllgridlha) vs the stored fonll-local/run/nlo_<q>/nl_<q>_.out: identical rows and the matched
   FONLL cross section hdmv + G (hdrs - hdml), G = pT^2 / (pT^2 + 25 m^2) (fragmfonll xsecmtch, c = 5), both as a
   relative difference and in units of FONLL's own integration error G * ehrs;
3. fragmfonll outputs vs the stored ones (bottom: central stream; charm: direct pseudoscalar and D* vector);
4. CTEQ6.6 inputs (step 02): charm calibration raw grid, fragmentation outputs and fitted weight; bottom
   public-point closure (raw rows vs the y = -0.5, 0, 0.5 rows of the stored 13-row grid; the point vs the public
   web form and the original local value).
Sections 2-4 and the stored-run parts compare with the original FONLL workspace of the optional reference tree
(GRENDEL_REF); without it they are skipped (a SKIP line). Section 1 and the CTEQ6.6 numbers always run.
"""
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np

REPRO = Path(os.environ["REPRO"])
sys.path.insert(0, str(REPRO / "pipelines/common"))
from ref import have_ref  # noqa: E402

WORK = Path(os.environ["FONLL_WORK"])
OUT = Path(os.environ["GRENDEL_FONLL_GRID_DIR"])
RUN = Path(os.environ["FONLL_RUN_DIR"])
TRACKED = Path(os.environ["FONLL_TRACKED_DIR"])
ORIG = Path(os.environ["FONLL_ORIG_WS"])
VAL = WORK / "validation"
VAL.mkdir(parents=True, exist_ok=True)

NAME = "fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_{}.dat"
MASS = {"bottom": 4.75, "charm": 1.50}
CPARAM = 5.0
# SHA-256 of the data rows (non-'#' lines) of every stored copy of the central grids (80 copies across the
# archive: fonll-local, fonll-nnpdf40 v0.1.0, the model repos, grendel-sensitivity at every tag)
DATA_SHA = {"bottom": "ec33ab8bc3e55ee02479a557a23ad2f11c6af857255850c5b77d65686f6536f2",
            "charm": "feb64ed416655d4a489c52cfb765e4ea41c62ce5ed3592605d8d46c89185bfb8"}
ORIG_WEIGHT = 1.2152533641182024

report: dict = {}
lines: list[str] = []


def say(msg: str = "") -> None:
    print(msg)
    lines.append(msg)


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def data_rows(path: Path) -> list[str]:
    return [l for l in path.read_text().splitlines() if not l.startswith("#")]


def rel_stats(a: np.ndarray, b: np.ndarray) -> dict:
    """|a - b| / |b| over b != 0; a and b must agree exactly where b == 0."""
    nz = b != 0
    zero_ok = bool(np.all(a[~nz] == 0))
    r = np.abs(a[nz] - b[nz]) / np.abs(b[nz])
    return {"n": int(nz.sum()), "max": float(r.max()) if r.size else 0.0,
            "p99": float(np.percentile(r, 99)) if r.size else 0.0,
            "median": float(np.median(r)) if r.size else 0.0, "zeros_match": zero_ok,
            "argmax": int(np.flatnonzero(nz)[r.argmax()]) if r.size else -1}


def sigma(pt, y, d2):
    return float(np.trapezoid(np.trapezoid(d2, y, axis=1), pt))


# ---------------------------------------------------------------- 1. meson grids
def meson(q: str, src: Path = OUT, tag: str = "") -> None:
    new, ref = src / NAME.format(q), TRACKED / NAME.format(q)
    key, quark = f"meson_{q}{tag}", q
    q = f"{q}{tag}"
    if not new.exists():
        say(f"[{q}] meson grid: not generated yet ({new})")
        return
    r: dict = {"file": str(new), "sha256": sha(new.read_bytes()), "tracked_sha256": sha(ref.read_bytes())}
    r["byte_identical"] = r["sha256"] == r["tracked_sha256"]
    hn = [l for l in new.read_text().splitlines() if l.startswith("#")]
    hr = [l for l in ref.read_text().splitlines() if l.startswith("#")]
    r["header_only_in_new"] = [l for l in hn if l not in hr]
    r["header_only_in_tracked"] = [l for l in hr if l not in hn]
    dn, dr = data_rows(new), data_rows(ref)
    r["data_sha256"] = sha(("\n".join(dn) + "\n").encode())
    r["data_identical"] = r["data_sha256"] == DATA_SHA[quark]
    r["identical_rows"] = sum(a == b for a, b in zip(dn, dr))
    r["rows"] = len(dr)
    a = np.loadtxt(new, comments="#")
    b = np.loadtxt(ref, comments="#")
    assert np.array_equal(a[:, :2], b[:, :2]), "node mismatch"
    st = rel_stats(a[:, 2], b[:, 2])
    i = st.pop("argmax")
    st["at_pT_y"] = [float(b[i, 0]), float(b[i, 1])] if i >= 0 else None
    r["per_bin_rel_diff"] = st
    pt, y = np.unique(b[:, 0]), np.unique(b[:, 1])
    an, bn = a[:, 2].reshape(len(pt), len(y)), b[:, 2].reshape(len(pt), len(y))
    s_new, s_ref = sigma(pt, y, an), sigma(pt, y, bn)
    r["sigma_pb"] = {"regenerated": s_new, "tracked": s_ref, "rel_diff": (s_new - s_ref) / s_ref}
    slices = {}
    for lo, hi in [(0, 5), (5, 10), (10, 20), (20, 50)]:
        m = (pt >= lo) & (pt <= hi)
        sn, sr = sigma(pt[m], y, an[m]), sigma(pt[m], y, bn[m])
        slices[f"pT {lo}-{hi}"] = (sn - sr) / sr
    for lo, hi in [(0, 1), (1, 2), (2, 3)]:
        m = (np.abs(y) >= lo) & (np.abs(y) <= hi)
        sn = sigma(pt, y[m], an[:, m]) if m.sum() > 1 else 0
        sr = sigma(pt, y[m], bn[:, m]) if m.sum() > 1 else 1
        slices[f"|y| {lo}-{hi} (one side)"] = (sn - sr) / sr if m.sum() > 1 else None
    r["sigma_slices_rel_diff"] = slices
    report[key] = r
    say(f"[{q}] meson grid vs tracked: byte-identical={r['byte_identical']}; data rows identical="
        f"{r['data_identical']} ({r['identical_rows']}/{r['rows']} rows equal)")
    say(f"[{q}]   header lines differing: {len(r['header_only_in_new'])} new / {len(r['header_only_in_tracked'])} tracked")
    say(f"[{q}]   per-bin |rel diff|: max {st['max']:.3e} at (pT, y) = {st['at_pT_y']}, p99 {st['p99']:.3e}, "
        f"median {st['median']:.3e}; pT=0 rows zero in both: {st['zeros_match']}")
    say(f"[{q}]   sigma (trapezoid, pT 0-50, |y|<3): {s_new:.9e} vs {s_ref:.9e} pb, rel {r['sigma_pb']['rel_diff']:+.3e}")
    say(f"[{q}]   slices: " + ", ".join(f"{k} {v:+.1e}" for k, v in slices.items() if v is not None))


# ---------------------------------------------------------------- 2. raw quark grids
def raw_cmp(label: str, new: Path, ref: Path, quark: str, rows=None) -> None:
    if not new.exists() or not ref.exists():
        say(f"[{label}] raw grid: missing ({new if not new.exists() else ref})")
        return
    tn, tr = new.read_text().splitlines(), ref.read_text().splitlines()
    npt = int(tr[0].split()[-1])
    body_r = tr[1:] if rows is None else [tr[1 + i * npt + j] for i in rows for j in range(npt)]
    body_n = tn[1:]
    if len(body_n) != len(body_r):
        say(f"[{label}] raw grid incomplete: {len(body_n)} of {len(body_r)} rows")
        return
    r = {"header_equal": tn[0] == tr[0], "identical_rows": sum(a == b for a, b in zip(body_n, body_r)),
         "rows": len(body_r)}
    a = np.array([[float(x) for x in l.split()] for l in body_n])
    b = np.array([[float(x) for x in l.split()] for l in body_r])
    assert np.array_equal(a[:, :2], b[:, :2])
    g = b[:, 0] ** 2 / (b[:, 0] ** 2 + CPARAM ** 2 * MASS[quark] ** 2)
    xa = a[:, 2] + g * (a[:, 4] - a[:, 3])
    xb = b[:, 2] + g * (b[:, 4] - b[:, 3])
    err = g * np.abs(b[:, 5])
    st = rel_stats(xa, xb)
    st.pop("argmax")
    r["matched_rel_diff"] = st
    pull = np.abs(xa - xb)[err > 0] / err[err > 0]
    r["matched_diff_over_fonll_error"] = {"max": float(pull.max()), "median": float(np.median(pull))}
    r["per_column_max_rel_diff"] = {name: rel_stats(a[:, c], b[:, c])["max"]
                                    for c, name in [(2, "hdmv"), (3, "hdml"), (4, "hdrs"), (5, "ehrs")]}
    report[f"raw_{label}"] = r
    say(f"[{label}] raw quark grid vs stored macOS: {r['identical_rows']}/{r['rows']} rows byte-identical; matched "
        f"FONLL xsec |rel diff| max {st['max']:.2e}, median {st['median']:.1e}; |diff|/(FONLL integration error) "
        f"max {r['matched_diff_over_fonll_error']['max']:.2e}")
    say(f"[{label}]   per column max rel diff: " +
        ", ".join(f"{k} {v:.1e}" for k, v in r["per_column_max_rel_diff"].items()))


def frag_cmp(label: str, new: Path, ref: Path) -> None:
    if not new.exists() or not ref.exists():
        say(f"[{label}] fragmfonll output: missing")
        return
    tn, tr = new.read_text().splitlines()[1:], ref.read_text().splitlines()[1:]
    r = {"identical_rows": sum(a == b for a, b in zip(tn, tr)), "rows": len(tr),
         "byte_identical": new.read_bytes() == ref.read_bytes()}
    if len(tn) == len(tr):
        a = np.array([[float(x.replace("D", "E")) for x in l.split()[:3]] for l in tn])
        b = np.array([[float(x.replace("D", "E")) for x in l.split()[:3]] for l in tr])
        st = rel_stats(a[:, 2], b[:, 2])
        st.pop("argmax")
        r["rel_diff"] = st
    report[f"frag_{label}"] = r
    msg = f"max |rel diff| {r['rel_diff']['max']:.2e}" if "rel_diff" in r else "row count differs"
    say(f"[{label}] fragmfonll output vs stored: {r['identical_rows']}/{r['rows']} rows identical; {msg}")


for q in ("bottom", "charm"):
    meson(q)
say()
HAVE_ORIG = have_ref("FONLL step 05: comparison with the stored original (macOS) FONLL run", ORIG / "run", ORIG / "output")
o = ORIG / "run"
if HAVE_ORIG:
    # xcheck_fragment.py: the stored macOS raw grids fragmented on this machine (isolates the steps after fonllgridlha)
    X = WORK / "xcheck"
    for q in ("bottom", "charm"):
        meson(q, X / "output", " xcheck(macOS raw grid, Linux fragmentation)")
    xr = X / "run"
    frag_cmp("xcheck bottom central", xr / "nlo_bottom_central" / "fragmfonll_central.dat",
             ORIG / "run" / "nlo_bottom" / "fragmfonll_central.dat")
    frag_cmp("xcheck charm direct", xr / "nlo_charm_central" / "fragmfonll_direct_pseudoscalar.dat",
             ORIG / "run" / "nlo_charm" / "fragmfonll_direct_pseudoscalar.dat")
    frag_cmp("xcheck charm D* vector", xr / "nlo_charm_central" / "fragmfonll_dstar_vector_feeddown.dat",
             ORIG / "run" / "nlo_charm" / "fragmfonll_dstar_vector_feeddown.dat")
    say()
    raw_cmp("bottom NNPDF4.0", RUN / "nlo_bottom_central" / "nl_b_.out", o / "nlo_bottom" / "nl_b_.out", "bottom")
    raw_cmp("charm NNPDF4.0", RUN / "nlo_charm_central" / "nl_c_.out", o / "nlo_charm" / "nl_c_.out", "charm")
    frag_cmp("bottom central", RUN / "nlo_bottom_central" / "fragmfonll_central.dat",
             o / "nlo_bottom" / "fragmfonll_central.dat")
    frag_cmp("charm direct", RUN / "nlo_charm_central" / "fragmfonll_direct_pseudoscalar.dat",
             o / "nlo_charm" / "fragmfonll_direct_pseudoscalar.dat")
    frag_cmp("charm D* vector", RUN / "nlo_charm_central" / "fragmfonll_dstar_vector_feeddown.dat",
             o / "nlo_charm" / "fragmfonll_dstar_vector_feeddown.dat")
    say()
    raw_cmp("charm CTEQ6.6 calibration", RUN / "validate_cteq66_charm_y0" / "val_c_.out",
            o / "validate_cteq66_charm_y0" / "val_c_.out", "charm")
    for s in ("direct_pseudoscalar", "dstar_vector"):
        frag_cmp(f"calibration {s}", RUN / "validate_cteq66_charm_y0" / f"fragmfonll_calibration_{s}.dat",
                 o / "validate_cteq66_charm_y0" / f"fragmfonll_calibration_{s}.dat")
else:
    say("SKIP the comparisons with the original (macOS) FONLL run, sections 2-4: no reference tree (GRENDEL_REF)")
    say()
cal = OUT / "charm_feeddown_calibration.json"
if cal.exists():
    c = json.loads(cal.read_text())
    w = c["vector_to_direct_weight"]
    report["charm_calibration"] = {"weight": w, "original_weight": ORIG_WEIGHT, "rel_diff": (w - ORIG_WEIGHT) / ORIG_WEIGHT,
                                   "max_abs_rel_vs_public": c["max_abs_relative_difference"]}
    say(f"[charm calibration] weight {w:.13g} vs original {ORIG_WEIGHT:.13g} (rel {(w - ORIG_WEIGHT) / ORIG_WEIGHT:+.2e}); "
        f"max |rel| vs the 9 public D0 points {c['max_abs_relative_difference']:.4e} (original 1.9410e-03)")
if HAVE_ORIG:
    raw_cmp("bottom CTEQ6.6 closure (y=-0.5,0,0.5)", RUN / "cteq66_bottom_y0" / "ct_b_.out",
            o / "cteq66_bottom" / "ct_b_.out", "bottom", rows=[5, 6, 7])
clo = VAL / "bottom_cteq66_public_point.json"
if clo.exists():
    c = json.loads(clo.read_text())
    if HAVE_ORIG:
        # the original point at full precision: fonll-local/run/cteq66_bottom/fragmfonll.dat, row (5, 0)
        for l in (ORIG / "run" / "cteq66_bottom" / "fragmfonll.dat").read_text().splitlines()[1:]:
            v = [float(x.replace("D", "E")) for x in l.split()]
            if v[0] == 5.0 and v[1] == 0.0:
                c["original_local_macos_full_precision"] = v[2]
                c["rel_vs_original_full_precision"] = (c["local_linux"] - v[2]) / v[2]
        orig = (f"original macOS {c.get('original_local_macos_full_precision', 0):.10e} "
                f"({c.get('rel_vs_original_full_precision', float('nan')):+.1e})")
    else:  # the 7-digit original value of data/fonll/PROVENANCE.md (closure.py)
        orig = f"original macOS {c['original_local_macos']:.6e} ({c['rel_vs_original_local']:+.1e})"
    report["bottom_public_point"] = c
    say(f"[bottom closure] B hadron pT=5 y=0 CTEQ6.6: {c['local_linux']:.7e} pb/GeV; public {c['public_fonll_web']:.5e} "
        f"({c['rel_vs_public']:+.3%}); {orig}")

(VAL / "compare.json").write_text(json.dumps(report, indent=2) + "\n")
(VAL / "summary.txt").write_text("\n".join(lines) + "\n")
print(f"\nwrote {VAL / 'compare.json'} and {VAL / 'summary.txt'}")
