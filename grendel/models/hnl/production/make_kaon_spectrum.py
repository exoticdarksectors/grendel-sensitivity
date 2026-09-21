#!/usr/bin/env python3
"""Rebuild the committed charged-kaon SoftQCD spectrum from tracked sources.

The sampler in ``generate_kaon_csvs.py`` reads a committed histogram,
``production/data/kaon_softqcd_spectrum.npz``. That histogram is produced by the
tracked Pythia driver ``kaon_softqcd.cc``; this script is the tracked recipe that
turns one into the other, so the committed ``.npz`` is reproducible rather than an
opaque binary.

Pipeline:
  1. compile ``kaon_softqcd.cc`` against Pythia 8.315 (via ``pythia8-config``) into
     a throwaway binary in a temp dir;
  2. run it for ``--n-events`` at ``--seed`` -- it streams final-state K+- ``(pT, y)``
     to stdout and run metadata (sigma_inel, <n_K+->) to stderr;
  3. bin ``(pT, y)`` into the committed ``200 x 160`` histogram (pT in ``[0, 10]`` GeV,
     y in ``[-8, 8]``) and write the ``.npz`` with the exact schema the sampler reads.

The committed spectrum was produced with the defaults below (Pythia 8.315,
n_events = 60000, seed = 42). Reproduce or check it with::

    python -m production.decay_engine.make_kaon_spectrum            # rebuild in place
    python -m production.decay_engine.make_kaon_spectrum --verify   # rebuild, compare, do NOT write

``--verify`` exits non-zero if the freshly generated histogram and metadata do not
match the committed file, so it doubles as a reproducibility regression check.
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from grendel.io.paths import third_party_dir

HERE = Path(__file__).resolve().parent
SRC = HERE / "kaon_softqcd.cc"
DATA = HERE / "data" / "kaon_softqcd_spectrum.npz"

N_EVENTS = 60_000
SEED = 42
PT_MAX, PT_BINS = 10.0, 200
Y_ABS, Y_BINS = 8.0, 160
# Fixed by kaon_softqcd.cc; embedded in the .npz for provenance and checked by --verify.
PROCESS = "SoftQCD:inelastic"
ECM_GEV = 14000.0
EXPECTED_VERSION = "8.315"   # the committed spectrum's Pythia; 8.317 shifts <n_K+-> ~0.5%


def pythia8_config():
    """Locate pythia8-config: ``$GRENDEL_PYTHIA8_CONFIG``, a Pythia 8.315 under
    the third-party directory, or one on PATH.

    The committed spectrum reproduces ONLY against 8.315: 8.317 (and
    presumably other minor versions) shifts <n_K+-> at the ~0.5% level, so a
    PATH fallback warns and ``--verify`` may then legitimately mismatch.
    """
    env = os.environ.get("GRENDEL_PYTHIA8_CONFIG")
    if env:
        return env
    cand = third_party_dir() / "pythia8315" / "bin" / "pythia8-config"
    if cand.exists():
        return str(cand)
    found = shutil.which("pythia8-config")
    if found:
        print(f"WARNING: Pythia 8.315 not found under {third_party_dir()}; using {found} -- "
              "this may not reproduce the committed 8.315 spectrum", file=sys.stderr)
        return found
    raise SystemExit("pythia8-config not found: set $GRENDEL_PYTHIA8_CONFIG or put one on PATH")
def _cfg(cfg, flag):
    return subprocess.check_output([cfg, flag], text=True).split()


def build(cfg, workdir):
    binary = workdir / "kaon_softqcd"
    cmd = ["c++", str(SRC), "-o", str(binary), *_cfg(cfg, "--cxxflags"), *_cfg(cfg, "--libs")]
    subprocess.run(cmd, check=True)
    return binary


def run(binary, cfg, n_events, seed):
    env = dict(os.environ)
    xmldoc = subprocess.check_output([cfg, "--xmldoc"], text=True).strip()
    env["PYTHIA8DATA"] = xmldoc  # locate the xmldoc at runtime regardless of cwd
    proc = subprocess.run(
        [str(binary), str(n_events), str(seed)],
        capture_output=True, text=True, env=env, check=True,
    )
    meta = {}
    for line in proc.stderr.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0].isupper():
            meta[parts[0]] = parts[1]
    pt, y = [], []
    for line in proc.stdout.splitlines():
        if not line.startswith("K "):
            continue
        _, a, b = line.split()
        pt.append(float(a))
        y.append(float(b))
    return np.asarray(pt), np.asarray(y), meta


def build_arrays(pt, y, meta, pythia_version):
    pt_edges = np.linspace(0.0, PT_MAX, PT_BINS + 1)
    y_edges = np.linspace(-Y_ABS, Y_ABS, Y_BINS + 1)
    hist, _, _ = np.histogram2d(pt, y, bins=[pt_edges, y_edges])
    n_kaon = int(meta["N_KAON"])
    n_events = int(meta["N_EVENTS"])
    if int(hist.sum()) != n_kaon:
        raise SystemExit(
            f"binned {int(hist.sum())} kaons but driver reported {n_kaon}; "
            "the (pT, y) window and the histogram edges must agree"
        )
    # Derive the normalization from the raw counts: the driver prints
    # NKAON_PER_EVT/SIGMA_KAON_PB to stderr at reduced precision, so recompute them
    # here to full float precision (this is how the committed spectrum was made).
    sigma_inel_mb = float(meta["SIGMA_INEL_MB"])
    n_kaon_per_inelastic = n_kaon / n_events
    sigma_kaon_pb = sigma_inel_mb * 1.0e9 * n_kaon_per_inelastic
    # Key order matches the committed .npz (savez_compressed writes in this order).
    return {
        "hist": hist,
        "pt_edges": pt_edges,
        "y_edges": y_edges,
        "sigma_inel_mb": np.float64(sigma_inel_mb),
        "n_kaon_per_inelastic": np.float64(n_kaon_per_inelastic),
        "sigma_kaon_pb": np.float64(sigma_kaon_pb),
        "n_events": np.int64(n_events),
        "n_kaon": np.int64(n_kaon),
        "pythia_version": np.asarray(pythia_version),
        "process": np.asarray(PROCESS),
        "ecm_gev": np.float64(ECM_GEV),
    }


def compare_to_committed(arrays):
    if not DATA.exists():
        raise SystemExit(f"no committed spectrum at {DATA} to verify against")
    ref = np.load(DATA)
    ok = True
    for key in ref.files:  # every committed field, incl. edges + version/process/ecm
        if key not in arrays:
            print(f"  {key} MISSING from rebuild", file=sys.stderr)
            ok = False
            continue
        if not np.array_equal(arrays[key], ref[key]):
            if key == "hist":
                detail = f"{int(np.count_nonzero(arrays[key] != ref[key]))}/{arrays[key].size} bins changed"
            else:
                detail = f"got {arrays[key]!r}, committed {ref[key]!r}"
            print(f"  {key} DIFFERS: {detail}", file=sys.stderr)
            ok = False
    extra = [k for k in arrays if k not in ref.files]
    if extra:
        print(f"  rebuild has EXTRA keys {extra}", file=sys.stderr)
        ok = False
    return ok


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n-events", type=int, default=N_EVENTS)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("-o", "--output", type=Path, default=DATA)
    ap.add_argument("--verify", action="store_true",
                    help="rebuild and compare to the committed .npz; do not write")
    ap.add_argument("--force", action="store_true",
                    help="allow writing even when the Pythia version is not "
                         f"{EXPECTED_VERSION} (which would not reproduce the committed spectrum)")
    args = ap.parse_args(argv)

    cfg = pythia8_config()
    version = subprocess.check_output([cfg, "--version"], text=True).strip()
    if not args.verify and version != EXPECTED_VERSION and not args.force:
        raise SystemExit(
            f"refusing to overwrite {args.output}: built against Pythia {version}, but the "
            f"committed spectrum is Pythia {EXPECTED_VERSION} (8.317 shifts <n_K+-> ~0.5%). "
            "Pass --verify to compare without writing, or --force to override."
        )
    with tempfile.TemporaryDirectory() as tmp:
        binary = build(cfg, Path(tmp))
        print(f"running {args.n_events} events (seed {args.seed}, Pythia {version}) ...",
              file=sys.stderr)
        pt, y, meta = run(binary, cfg, args.n_events, args.seed)
    arrays = build_arrays(pt, y, meta, version)
    print(
        f"sigma_inel={float(arrays['sigma_inel_mb']):.4f} mb, "
        f"<n_K+->={float(arrays['n_kaon_per_inelastic']):.4f}, "
        f"n_kaon={int(arrays['n_kaon'])}, "
        f"sigma_kaon={float(arrays['sigma_kaon_pb']):.4e} pb",
        file=sys.stderr,
    )

    if args.verify:
        ok = compare_to_committed(arrays)
        print("VERIFY: match" if ok else "VERIFY: MISMATCH", file=sys.stderr)
        return 0 if ok else 1

    np.savez_compressed(args.output, **arrays)
    print(f"wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
