#!/usr/bin/env python3
"""BC5 decay templates with the Winkler lifetime and channel split and exHad final states.

Each exHad bundle (``exhad:scalar-1809``, the v2 scalar templates) is copied with
- ``ctau_m_u2eq1`` set to the Winkler model value the bundle stores as
  ``ctau_m_u2eq1_reference`` (equal to ``model.ctau_sin2theta1(m, 'winkler')``);
- its events subsampled, without replacement, so that the fractions of the
  channel groups e e, mu mu, tau tau and hadronic are the Winkler branching
  ratios of ``model.branching_ratios(m, 'winkler')`` (hadronic = every
  non-leptonic channel). Inside the hadronic group the composition and every
  final state stay exHad's.
The acceptance MC draws templates uniformly, so the subsampled bundle carries
the Winkler split directly.

    python make_hybrid_templates.py SRC_DIR DST_DIR
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from grendel.models.scalar import model  # noqa: E402

GROUP_OF_LABEL = {"ePeM": "ee", "muPmuM": "mumu", "tauPtauM": "tautau", "hadronic": "hadronic"}
LEPTONIC = ("ee", "mumu", "tautau")
SEED = 20261007
PER_DAUGHTER = ("pdg", "px", "py", "pz", "energy", "mass", "charge", "stable")
PER_EVENT = ("daughter_counts", "channel_label", "n_charged")


def winkler_groups(m):
    br = model.branching_ratios(m, scheme="winkler")
    out = {g: float(br.get(g, 0.0)) for g in LEPTONIC}
    out["hadronic"] = float(sum(v for k, v in br.items() if k not in LEPTONIC))
    total = sum(out.values())
    return {g: v / total for g, v in out.items()}


def subsample(d, m, rng):
    groups = np.array([GROUP_OF_LABEL[str(x)] for x in d["channel_label"]])
    target = winkler_groups(m)
    have = {g: int((groups == g).sum()) for g in target}
    missing = {g: t for g, t in target.items() if t > 0 and have[g] == 0}
    if any(t > 1e-3 for t in missing.values()):
        raise SystemExit(f"m={m}: exHad bundle has no events for {missing}")
    target = {g: (0.0 if g in missing else t) for g, t in target.items()}
    norm = sum(target.values())
    target = {g: t / norm for g, t in target.items()}
    n_total = min(have[g] / t for g, t in target.items() if t > 0)
    keep = []
    for g, t in target.items():
        n_g = int(round(t * n_total))
        idx = np.flatnonzero(groups == g)
        if n_g:
            keep.append(rng.choice(idx, size=min(n_g, len(idx)), replace=False))
    keep = np.sort(np.concatenate(keep))
    counts = np.asarray(d["daughter_counts"])
    off = np.concatenate([[0], np.cumsum(counts)])
    flat = np.concatenate([np.arange(off[i], off[i + 1]) for i in keep])
    out = dict(d)
    for k in PER_EVENT:
        out[k] = np.asarray(d[k])[keep]
    for k in PER_DAUGHTER:
        out[k] = np.asarray(d[k])[flat]
    out["n_templates"] = np.int32(len(keep))
    got = {g: float((groups[keep] == g).mean()) for g in target}
    return out, target, got, have


def main(src, dst):
    src, dst = Path(src), Path(dst)
    dst.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    records = []
    for f in sorted(src.glob("templates_*.npz")):
        d = dict(np.load(f, allow_pickle=False))
        m = float(d["mass_GeV"])
        ref = float(d["ctau_m_u2eq1_reference"])
        winkler = float(model.ctau_sin2theta1(m, scheme="winkler"))
        if not np.isfinite(ref) or ref <= 0 or abs(ref / winkler - 1) > 1e-6:
            raise SystemExit(f"{f.name}: reference ctau {ref} is not the Winkler model value {winkler}")
        out, target, got, have = subsample(d, m, rng)
        out["ctau_m_u2eq1_exhad"] = d["ctau_m_u2eq1"]
        out["ctau_m_u2eq1"] = np.float64(ref)
        out["ctau_source"] = np.array("Winkler model (ctau_m_u2eq1_reference of the exHad bundle)")
        out["decay_model"] = np.array(str(d["decay_model"]) + "+winkler-ctau-br")
        np.savez_compressed(dst / f.name, **out)
        records.append({"file": f.name, "mass_GeV": m, "n_templates": int(out["n_templates"]),
                        "ctau_m_u2eq1": ref, "ctau_m_u2eq1_exhad": float(d["ctau_m_u2eq1"]),
                        "exhad_events": have, "winkler_fractions": target, "kept_fractions": got})
        print(f"{f.name}: {int(out['n_templates'])} templates, "
              + ", ".join(f"{g} {got[g]:.3f}/{target[g]:.3f}" for g in target if target[g] > 0))
    (dst / "MANIFEST.json").write_text(json.dumps(
        {"variant": "exHad final states (scalar-1809) with the Winkler lifetime and e e / mu mu / "
                    "tau tau / hadronic split", "source_templates": str(src), "seed": SEED,
         "points": records}, indent=1))
    print(f"wrote {len(records)} bundles -> {dst}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
