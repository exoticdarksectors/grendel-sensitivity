"""Generate exHad rest-frame decay templates and cache them as ``.npz`` bundles.

One driver for every model. A model supplies a ``TemplateModel``: which exHad
model to run, the points to generate (mass, and flavor/mixing for the HNL),
where each bundle goes, how the lifetime at unit coupling is read off exHad's
tables, its own reference lifetime for the diagnostics, any extra bundle
fields, and the seed policy. The driver owns the exHad session, the retry
strategy, the bundle layout the acceptance Monte Carlo reads, and the
manifest.

Bundle layout (the analysis reads only the first group)::

    templates_<label>.npz
        daughter_counts (N,)            int32   daughters per template
        pdg             (M,)            int32   M = sum(daughter_counts)
        px, py, pz, energy, mass (M,)   float64 rest-frame 4-momentum + mass
        charge          (M,)            float64 in units of e
        stable          (M,)            bool    always True (exHad returns final particles)
        mass_GeV, ctau_m_u2eq1, n_templates, seed   scalars
        flavor                          str
      provenance / diagnostics:
        channel_label (N,), n_charged (N,), decay_model, exhad_commit,
        exhad_variation, pythia_version, exhad_seed, ctau_source,
        ctau_m_u2eq1_reference, and the model's extra fields

exHad leaves mu, pi+-, K+-, K0_L and n undecayed and decays everything else
(pi0, K0_S, Lambda, tau, ...) at the vertex. The few-cm flight of a K0_S or a
hyperon is negligible for a detector 22 m away.

Run under exHad's own virtual environment (it has the ``exhad`` package,
numpy and ``particle``); ``--exhad-root`` defaults to the checkout made by
``third_party/fetch.py``.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import platform
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np

from ..io.paths import third_party_dir
from .templates import (channel_fractions, complete_manifest_points, flatten_events,
                        generate_all_robust, generation_runs, import_exhad,
                        load_prior_manifest, point_seed, read_exhad_provenance)


def default_exhad_root() -> Path:
    return third_party_dir() / "exHad"


@dataclass(frozen=True)
class TemplatePoint:
    mass: float
    label: str
    flavor: str | None = None          # HNL only
    mixing: tuple | None = None        # HNL only: exHad mixing vector

    @property
    def tag(self) -> str:
        return f"{self.flavor} m={self.mass:.3f}" if self.flavor else f"m={self.mass:.3f}"


@dataclass
class TemplateModel:
    """What a model tells the generator."""
    exhad_model: str                       # exHad model name (``hnl``, ``scalar-1809``, ``alp-fermion``)
    decay_model: str                       # written to the bundle and manifest
    flavor_tag: str | None                 # bundle ``flavor`` for flavorless models (``BC4``, ``BC10``)
    seed_namespace: str
    ctau_source: str
    ctau_ref: Callable[[TemplatePoint], float]          # exHad lifetime at unit coupling [m]
    reference_ctau: Callable[[TemplatePoint], float]    # the model's own lifetime (diagnostic), NaN if none
    dest: Callable[[TemplatePoint, Path], Path]
    strategy: str = "sub_request"          # "sub_request" (robust, seeded sub-requests) or "attempt" (retry whole point with seed+attempt)
    extra_bundle: Callable[[TemplatePoint, dict], dict] = lambda pt, bundle: {}
    reference_source: str | None = None
    manifest_extra: dict = field(default_factory=dict)
    skip: Callable[[TemplatePoint], str | None] = lambda pt: None   # reason to skip a point

    def seed(self, base_seed: int, pt: TemplatePoint) -> int:
        parts = [pt.flavor, pt.label] if pt.flavor else [pt.label]
        if self.exhad_model.startswith("scalar"):
            parts = [self.exhad_model, pt.label]
        return point_seed(base_seed, *parts, namespace=self.seed_namespace)


def add_generation_options(ap: argparse.ArgumentParser, *, robust: bool) -> None:
    ap.add_argument("--n-templates", type=int, default=20_000)
    ap.add_argument("--out", type=Path, default=None,
                    help="template directory (default: the model's templates path)")
    ap.add_argument("--seed", type=int, default=1234,
                    help="base seed; every point gets a distinct seed derived from it")
    ap.add_argument("--workers", type=int, default=None,
                    help="exHad worker processes (default: exHad's own)")
    ap.add_argument("--chunk-size", type=int, default=512)
    if robust:
        ap.add_argument("--sub-request", type=int, default=2000,
                        help="events per exHad request; a failed request is retried alone with a fresh seed")
    ap.add_argument("--max-retries", type=int, default=3,
                    help="retries per failed request, each with a fresh derived seed")
    ap.add_argument("--variation", default="central", help="exHad model variation")
    ap.add_argument("--exhad-root", type=Path, default=default_exhad_root())
    ap.add_argument("--skip-existing", action="store_true")


def _generate(model, generator, pt, n_templates, seed, args):
    """(events, labels, retries, seed_used) with the model's retry strategy."""
    if model.strategy == "attempt":
        last_error = None
        for attempt in range(args.max_retries):
            try:
                kwargs = {"mixing": list(pt.mixing)} if pt.mixing is not None else {}
                sample = generator.generate_all(float(pt.mass), int(n_templates),
                                                seed=seed + attempt, **kwargs)
                events, labels = sample["events"], sample["channel_labels"]
                if len(events) != n_templates:
                    raise RuntimeError(f"exHad returned {len(events)} events for {n_templates} requested")
                return events, labels, [f"attempt {a + 1}" for a in range(attempt)], seed + attempt
            except RuntimeError as err:   # e.g. a narrow interval near a quark threshold
                last_error = err
                print(f"  [retry {attempt + 1}/{args.max_retries}] {pt.tag}: {str(err)[:160]}", flush=True)
        raise RuntimeError(str(last_error))
    events, labels, retries = generate_all_robust(
        generator, pt.mass, n_templates, seed,
        sub_request=args.sub_request, max_retries=args.max_retries)
    return events, labels, retries, seed


def generate_one(model, generator, pt, out_dir, base_seed, provenance, args):
    reason = model.skip(pt)
    if reason:
        print(f"  {pt.tag}: skip ({reason})", flush=True)
        return None
    dest = model.dest(pt, out_dir)
    if args.skip_existing and dest.exists():
        print(f"  [skip] {pt.tag} ({dest.name})", flush=True)
        return None
    seed = model.seed(base_seed, pt)
    t0 = time.time()
    try:
        events, labels, retries, seed_used = _generate(model, generator, pt, args.n_templates, seed, args)
    except RuntimeError as err:
        print(f"  [FAIL] {pt.tag}: {err}", flush=True)
        return {"flavor": pt.flavor, "mass_GeV": float(pt.mass), "error": str(err)}
    bundle = flatten_events(events, labels)
    ctau_m = float(model.ctau_ref(pt))
    ctau_ref = float(model.reference_ctau(pt))
    vis_frac = float((bundle["n_charged"] >= 2).mean())
    dest.parent.mkdir(parents=True, exist_ok=True)
    fields = dict(
        mass_GeV=np.float64(pt.mass),
        ctau_m_u2eq1=np.float64(ctau_m),
        n_templates=np.int32(args.n_templates),
        seed=np.int64(base_seed),
        flavor=np.array(pt.flavor or model.flavor_tag),
        decay_model=np.array(model.decay_model),
        exhad_commit=np.array(provenance["exhad_commit"]),
        exhad_variation=np.array(args.variation),
        pythia_version=np.array(provenance["pythia_version"]),
        exhad_seed=np.uint64(seed_used),
        ctau_source=np.array(model.ctau_source),
        ctau_m_u2eq1_reference=np.float64(ctau_ref),
    )
    if pt.mixing is not None:
        fields["mixing"] = np.asarray(pt.mixing, dtype=np.float64)
    if model.reference_source:
        fields["ctau_reference_source"] = np.array(model.reference_source)
    fields.update(model.extra_bundle(pt, bundle))
    np.savez_compressed(dest, **fields, **bundle)

    frac = channel_fractions(labels)
    top = ", ".join(f"{k}:{v:.0%}" for k, v in list(frac.items())[:3])
    ratio = f" (x{ctau_m / ctau_ref:.3f} vs reference)" if np.isfinite(ctau_ref) and ctau_ref > 0 else ""
    note = f", retries {len(retries)}" if retries else ""
    print(f"  [ok]   {pt.tag}  vis(>=2 chg)={vis_frac:5.1%}  <mult>={bundle['daughter_counts'].mean():.2f}  "
          f"ctau(ref)={ctau_m:.3e} m{ratio}  top: {top}  [{time.time() - t0:.1f}s{note}] -> {dest.name}",
          flush=True)
    record = {"mass_GeV": float(pt.mass), "file": str(dest.relative_to(out_dir)),
              "exhad_seed": int(seed_used), "retries": retries,
              "vis_frac_ge2_charged": vis_frac,
              "mean_multiplicity": float(bundle["daughter_counts"].mean()),
              "ctau_m_u2eq1": ctau_m,
              "ctau_m_u2eq1_reference": ctau_ref if np.isfinite(ctau_ref) else None,
              "channel_fractions": frac}
    if pt.flavor:
        record = {"flavor": pt.flavor, **record}
    return record


def run(model: TemplateModel, points: list[TemplatePoint], out_dir: Path, args) -> int:
    """Generate every point, then write ``MANIFEST.json``. Returns the exit code."""
    exhad_root = args.exhad_root.expanduser().resolve()
    Generator = import_exhad(exhad_root).Generator
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    provenance = read_exhad_provenance(exhad_root)
    print(f"Generating exHad {model.exhad_model} decay templates -> {out_dir}")
    print(f"  exHad {provenance['exhad_commit'][:12]} (Pythia {provenance['pythia_version']}), "
          f"variation={args.variation}, points={len(points)}, n_templates={args.n_templates}, "
          f"base seed={args.seed}", flush=True)

    # A resumed run regenerates only the missing points; the manifest must
    # still list every bundle, so the prior manifest's records are carried over.
    prior_manifest = load_prior_manifest(out_dir)
    records, failures = [], []
    t_start = time.time()
    with Generator(model.exhad_model, workers=args.workers, chunk_size=args.chunk_size,
                   variation=args.variation, root=str(exhad_root)) as generator:
        for pt in points:
            rec = generate_one(model, generator, pt, out_dir, args.seed, provenance, args)
            if rec is None:
                continue
            (failures if "error" in rec else records).append(rec)

    points_done, points_provenance = complete_manifest_points(
        out_dir, records, prior_manifest.get("points", []))
    this_run = {
        "generated_at": _dt.datetime.now().isoformat(timespec="seconds"),
        "n_generated": len(records), "n_failures": len(failures),
        "skip_existing": bool(args.skip_existing),
        "elapsed_s": round(time.time() - t_start, 1),
    }
    if model.strategy == "sub_request":
        this_run.update(sub_request=args.sub_request, max_retries=args.max_retries)
        seed_policy = (f"sha256('{model.seed_namespace}/<base>/<point>')[:8] mod 2^63; sub-requests of "
                       f"{args.sub_request} events seeded by sha256('exhad-sub-request/v1/<seed>/<index>/<attempt>')")
    else:
        seed_policy = f"sha256('{model.seed_namespace}/<base>/<point>')[:8] mod 2^63 (+attempt on retry)"
    manifest = {
        "decay_model": model.decay_model,
        "exhad_root": str(exhad_root),
        **provenance,
        "exhad_variation": args.variation,
        "ctau_source": model.ctau_source,
        **({"ctau_reference_source": model.reference_source} if model.reference_source else {}),
        "n_templates": args.n_templates,
        "base_seed": args.seed,
        "seed_policy": seed_policy,
        "chunk_size": args.chunk_size,
        **({"sub_request": args.sub_request} if model.strategy == "sub_request" else {}),
        "max_retries": args.max_retries,
        "n_points": len(points),
        **model.manifest_extra,
        "generated_at": this_run["generated_at"],
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "platform": platform.platform(),
        "elapsed_s": this_run["elapsed_s"],
        "points": points_done,
        "points_provenance": points_provenance,
        "generation_runs": generation_runs(prior_manifest, this_run),
        "failures": failures,
        "n_retries": int(sum(len(r.get("retries", [])) for r in points_done)),
    }
    manifest_path = out_dir / "MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"Wrote {len(records)} template files, {len(failures)} failures; manifest lists "
          f"{len(points_done)} bundles ({points_provenance}), {manifest['n_retries']} retries "
          f"-> {manifest_path}")
    return 1 if failures else 0
