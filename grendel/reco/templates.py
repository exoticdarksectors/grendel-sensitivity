"""Helpers for exHad-generated rest-frame decay templates.

Used by the per-model template generators. exHad (Kryshtal & Ovchynnikov,
arXiv:2609.16104) is fetched by ``third_party/fetch.py``; these helpers only
pack its events into the flat ``.npz`` bundle the acceptance Monte Carlo
(``grendel.reco.acceptance.build_event_mc``) reads, derive seeds, and read
provenance. They import nothing from exHad themselves, so they are testable
without its runtime.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np

HBARC_GEV_M = 1.973269804e-16  # hbar*c in GeV*m


def import_exhad(vendored_root):
    """Import the installed ``exhad`` package from exHad's virtual environment.

    ``vendored_root`` is the exHad checkout (``third_party/exHad``). If its
    parent directory is on ``sys.path`` a bare ``exhad`` directory there would
    be found first and taken as an empty namespace package, and ``from exhad
    import Generator`` would fail with "unknown location"; strip such entries
    for the duration of the import.
    """
    import importlib
    import sys
    parent = Path(vendored_root).resolve().parent
    kept = list(sys.path)
    sys.path[:] = [p for p in sys.path if Path(p or ".").resolve() != parent]
    try:
        module = sys.modules.get("exhad")
        if module is not None and getattr(module, "__file__", None) is None:
            del sys.modules["exhad"]  # a namespace-package import already happened
        exhad = importlib.import_module("exhad")
    finally:
        sys.path[:] = kept
    if getattr(exhad, "__file__", None) is None:
        raise ImportError(f"exhad resolved to a namespace package; run with {vendored_root}/.venv/bin/python")
    return exhad

_CHARGE_CACHE: dict[int, float] = {}


def charge_of(pdg: int) -> float:
    """Electric charge in units of e from the PDG code (via the ``particle`` package)."""
    pdg = int(pdg)
    if pdg not in _CHARGE_CACHE:
        from particle import Particle
        try:
            q = float(Particle.from_pdgid(pdg).charge)
        except Exception:  # unknown code: treat as neutral, but say so
            print(f"  WARNING: unknown PDG code {pdg}; charge set to 0", flush=True)
            q = 0.0
        _CHARGE_CACHE[pdg] = q
    return _CHARGE_CACHE[pdg]


def read_exhad_provenance(exhad_root: Path) -> dict:
    """Upstream commit (UPSTREAM_COMMIT.txt written by third_party/fetch.py, or the
    checkout's git HEAD) and the Pythia version the build used."""
    exhad_root = Path(exhad_root)
    commit_file = exhad_root / "UPSTREAM_COMMIT.txt"
    commit = "unknown"
    if commit_file.exists():
        commit = commit_file.read_text().strip()
    elif (exhad_root / ".git").exists():
        import subprocess
        try:
            commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=exhad_root, check=True,
                                    capture_output=True, text=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            pass
    pythia_version = "unknown"
    build = exhad_root / "cpp" / ".pythia8-build"
    if build.exists():
        for line in build.read_text().splitlines():
            if line.startswith("version="):
                pythia_version = line.split("=", 1)[1].strip()
    return {"exhad_commit": commit, "pythia_version": pythia_version}


def point_seed(base_seed: int, *parts, namespace: str) -> int:
    """Deterministic, distinct exHad seed per point; unsigned 63-bit.

    The key is ``<namespace>/<base_seed>/<part>/<part>...``; the HNL generator uses
    ``namespace="exhad-hnl-templates/v1"`` with parts ``(flavor, label)``.
    """
    key = "/".join([namespace, str(base_seed), *[str(p) for p in parts]]).encode()
    return int.from_bytes(hashlib.sha256(key).digest()[:8], "little") % (2 ** 63)


def flatten_events(events, labels):
    """Pack exHad events (``[[px,py,pz,E,m,pdg], ...]`` per decay) into the flat bundle.

    All returned particles are final (exHad decays everything but mu, pi+-, K+-, K0_L, n),
    so ``stable`` is all True. ``n_charged`` counts |charge| > 0.5 particles per decay.
    """
    counts = np.array([len(ev) for ev in events], dtype=np.int32)
    total = int(counts.sum())
    pdg = np.empty(total, dtype=np.int32)
    px = np.empty(total); py = np.empty(total); pz = np.empty(total)
    energy = np.empty(total); mass = np.empty(total); charge = np.empty(total)
    n_charged = np.zeros(len(events), dtype=np.int16)
    k = 0
    for i, ev in enumerate(events):
        nq = 0
        for part in ev:
            px[k], py[k], pz[k], energy[k], mass[k] = part[0], part[1], part[2], part[3], part[4]
            code = int(part[5])
            pdg[k] = code
            q = charge_of(code)
            charge[k] = q
            if abs(q) > 0.5:
                nq += 1
            k += 1
        n_charged[i] = nq
    return dict(
        daughter_counts=counts, pdg=pdg, px=px, py=py, pz=pz, energy=energy,
        mass=mass, charge=charge, stable=np.ones(total, dtype=bool),
        channel_label=np.array(list(labels), dtype=str), n_charged=n_charged,
    )


def sub_request_seed(seed: int, index: int, attempt: int) -> int:
    """Seed of sub-request ``index`` (retry ``attempt``) of a request seeded with ``seed``."""
    key = f"exhad-sub-request/v1/{seed}/{index}/{attempt}".encode()
    return int.from_bytes(hashlib.sha256(key).digest()[:8], "little") % (2 ** 63)


def generate_all_robust(generator, mass, n_events, seed, *, sub_request=2000,
                        max_retries=3, log=print, **kwargs):
    """``generator.generate_all`` in independent sub-requests with per-sub-request retries.

    exHad fails a whole request when one event fails a consistency check (observed as
    "Full decay does not conserve four-momentum", well below 1e-4 per event). Splitting a
    request keeps such a failure from discarding the point: only the failed sub-request is
    regenerated, with a fresh derived seed. The events of a sub-request depend on its own seed
    and size only, so the sample is a union of independent draws from the same distribution.

    Returns ``(events, channel_labels, retries)`` where ``retries`` lists the failed
    ``(index, attempt, message)`` triples.
    """
    events, labels, retries = [], [], []
    n_sub = max(1, math.ceil(n_events / sub_request))
    for index in range(n_sub):
        n_i = min(sub_request, n_events - index * sub_request)
        if n_i <= 0:
            break
        last = None
        for attempt in range(max_retries + 1):
            s = sub_request_seed(seed, index, attempt)
            try:
                out = generator.generate_all(float(mass), int(n_i), seed=s, **kwargs)
                break
            except RuntimeError as err:
                last = err
                retries.append((index, attempt, str(err)[:200]))
                log(f"  [retry] mass={mass} sub-request {index} attempt {attempt}: {str(err)[:120]}")
        else:
            raise RuntimeError(f"exHad sub-request {index} failed {max_retries + 1} times: {last}")
        events.extend(out["events"])
        labels.extend(out["channel_labels"])
    if len(events) != n_events:
        raise RuntimeError(f"exHad returned {len(events)} events for {n_events} requested")
    return events, labels, retries


def loglog_interp(x: float, xs, ys) -> float:
    """Interpolate a positive table linearly in log-log; raises outside the table."""
    xs = np.asarray(xs, float); ys = np.asarray(ys, float)
    if not (xs[0] <= x <= xs[-1]):
        raise ValueError(f"{x} outside the table [{xs[0]}, {xs[-1]}]")
    return float(np.exp(np.interp(np.log(x), np.log(xs), np.log(ys))))


def channel_fractions(labels) -> dict:
    counts: dict[str, int] = {}
    for lab in labels:
        counts[lab] = counts.get(lab, 0) + 1
    n = max(1, len(labels))
    return {k: v / n for k, v in sorted(counts.items(), key=lambda kv: -kv[1])}


# ---- manifest bookkeeping shared by the three generators ---------------------

def load_prior_manifest(out_dir) -> dict:
    """The MANIFEST.json a previous run left in ``out_dir``, or {} when there is none."""
    path = Path(out_dir) / "MANIFEST.json"
    if not path.exists():
        return {}
    import json
    return json.loads(path.read_text())


def record_from_bundle(path, out_dir) -> dict:
    """Manifest record of an existing template bundle, rebuilt from its own keys.

    Everything a record holds is stored in the bundle (seed, lifetime and its reference,
    channel labels, charged multiplicities) except the retry log of the request that
    produced it, so a rebuilt record says so in ``record_source``.
    """
    path = Path(path)
    with np.load(path, allow_pickle=False) as d:
        n_charged = np.asarray(d["n_charged"])
        counts = np.asarray(d["daughter_counts"])
        labels = d["channel_label"].tolist()
        reference = (float(d["ctau_m_u2eq1_reference"])
                     if "ctau_m_u2eq1_reference" in d.files else float("nan"))
        record = {
            "mass_GeV": float(d["mass_GeV"]),
            "file": str(path.relative_to(Path(out_dir))),
            "exhad_seed": int(d["exhad_seed"]),
            "vis_frac_ge2_charged": float((n_charged >= 2).mean()),
            "mean_multiplicity": float(counts.mean()),
            "ctau_m_u2eq1": float(d["ctau_m_u2eq1"]),
            "ctau_m_u2eq1_reference": reference if math.isfinite(reference) else None,
            "channel_fractions": channel_fractions(labels),
            "record_source": "rebuilt from the bundle; the retry log of its request was not retained",
        }
        flavor = str(d["flavor"]) if "flavor" in d.files else ""
    if flavor in ("Ue", "Umu", "Utau"):
        record = {"flavor": flavor, **record}
    return record


def complete_manifest_points(out_dir, generated, prior_points):
    """Manifest records for EVERY bundle under ``out_dir``.

    A resumed run (``--skip-existing``) generates only the missing points, so its own
    records cover a subset. The rest come from the prior manifest when it has a record
    for that bundle (it carries the retry log), otherwise from the bundle itself. Returns
    ``(points, counts)`` with the points sorted by (flavor, mass) and ``counts`` saying
    where each came from.
    """
    out_dir = Path(out_dir)
    this_run = {rec["file"] for rec in generated}
    prior = {rec["file"]: rec for rec in prior_points if "file" in rec and "error" not in rec}
    points = list(generated)
    counts = {"generated": len(generated), "from_prior_manifest": 0, "rebuilt_from_bundle": 0}
    for path in sorted(out_dir.rglob("templates_*.npz")):
        rel = str(path.relative_to(out_dir))
        if rel in this_run:
            continue
        if rel in prior:
            points.append(prior[rel])
            counts["from_prior_manifest"] += 1
        else:
            points.append(record_from_bundle(path, out_dir))
            counts["rebuilt_from_bundle"] += 1
    points.sort(key=lambda rec: (rec.get("flavor", ""), rec["mass_GeV"]))
    return points, counts


def generation_runs(prior_manifest: dict, this_run: dict) -> list:
    """The ``generation_runs`` history: the prior manifest's list (or, for a manifest
    written before this key existed, one entry reconstructed from its top-level fields),
    followed by ``this_run``."""
    runs = list(prior_manifest.get("generation_runs", []))
    if prior_manifest and not runs:
        import re
        policy = str(prior_manifest.get("seed_policy", ""))
        found = re.search(r"sub-requests of (\d+) events", policy)
        runs.append({
            "generated_at": prior_manifest.get("generated_at"),
            "n_generated": len(prior_manifest.get("points", [])),
            "sub_request": int(found.group(1)) if found else None,
            "max_retries": prior_manifest.get("max_retries"),
            "skip_existing": None,
            "elapsed_s": prior_manifest.get("elapsed_s"),
            "note": "reconstructed from a manifest written before generation_runs was recorded",
        })
    runs.append(this_run)
    return runs
