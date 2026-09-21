"""Legacy two-track proxy templates for the BC10 acceptance Monte Carlo.

Produces, per ALP mass, the same flat npz bundle the other template
generators write, so ``grendel.reco.acceptance.build_event_mc`` consumes it
unchanged:

    daughter_counts (N,) int32      pdg (M,) int32
    px, py, pz, energy, mass (M,)   charge (M,)   stable (M,) bool
    mass_GeV, ctau_m_u2eq1, n_templates, seed, flavor   scalars

Each template is one isotropic two-body rest-frame decay a -> X+ X- of a
channel drawn with probability BR(channel) (``model.visible_channel_weights``),
so uniform template sampling reproduces the visible-BR mix. The two charged
daughters are the two reconstructed tracks; for the multi-prong modes
(tau tau, hadronic) they proxy the two leading charged tracks. ``ctau_m_u2eq1``
is the lifetime at the reference coupling ``model.INV_F_REF``; the scan
rescales it by 1/u2.

This approximation is kept only to reproduce the earliest BC10 scan; the
published templates come from ``templates_pythia`` (full branching) or
``templates_exhad``.

    python -m grendel.models.alp.templates_proxy --legacy-proxy --n-templates 20000
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from ...io.paths import ModelPaths
from ...io.vectors import format_mass_for_filename
from . import model
from .mass_grid import ALP_MASS_GRID


def _isotropic_unit(n, rng):
    cos_t = rng.uniform(-1.0, 1.0, n)
    sin_t = np.sqrt(np.maximum(1.0 - cos_t ** 2, 0.0))
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    return np.column_stack([sin_t * np.cos(phi), sin_t * np.sin(phi), cos_t])


def build_templates(m_a, n_templates, rng):
    """Return the flat npz-ready dict for one ALP mass (or None if no open
    visible channel)."""
    weights = model.visible_channel_weights(m_a, model.INV_F_REF)
    if not weights:
        return None
    channels = list(weights)
    probs = np.array([weights[c] for c in channels], dtype=float)
    probs /= probs.sum()
    pick = rng.choice(len(channels), size=n_templates, p=probs)

    dirs = _isotropic_unit(n_templates, rng)
    pdg = np.empty(2 * n_templates, dtype=np.int32)
    px = np.empty(2 * n_templates); py = np.empty(2 * n_templates)
    pz = np.empty(2 * n_templates); energy = np.empty(2 * n_templates)
    mass = np.empty(2 * n_templates); charge = np.empty(2 * n_templates)

    for k in range(n_templates):
        ch = channels[pick[k]]
        pid, m_d = model.VISIBLE_CHANNELS[ch]
        if m_a <= 2.0 * m_d:                  # numerical guard at threshold
            m_d = 0.499 * m_a
        p_star = np.sqrt(max((m_a / 2.0) ** 2 - m_d ** 2, 0.0))
        e_star = m_a / 2.0
        d = dirs[k]
        i = 2 * k
        # daughter +  (pdg +pid, charge +1)
        px[i], py[i], pz[i] = p_star * d
        energy[i] = e_star; mass[i] = m_d; charge[i] = +1.0; pdg[i] = pid
        # daughter -  (pdg -pid, charge -1), back-to-back
        px[i + 1], py[i + 1], pz[i + 1] = -p_star * d
        energy[i + 1] = e_star; mass[i + 1] = m_d
        charge[i + 1] = -1.0; pdg[i + 1] = -pid

    return dict(
        daughter_counts=np.full(n_templates, 2, dtype=np.int32),
        pdg=pdg, px=px, py=py, pz=pz, energy=energy, mass=mass,
        charge=charge, stable=np.ones(2 * n_templates, dtype=bool),
        mass_GeV=np.float64(m_a),
        ctau_m_u2eq1=np.float64(model.alp_ctau(m_a, model.INV_F_REF)),
        n_templates=np.int32(n_templates),
        seed=np.int64(-1), flavor=np.array("BC10"),
    )


def generate(masses, n_templates, out_dir, seed=1234):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for m_a in masses:
        resonance = model.excluded_light_meson_resonance(m_a)
        if resonance is not None:
            print(f"  m_a={m_a:.3f}: skip unsupported {resonance} resonance")
            continue
        rng = np.random.default_rng(seed + int(round(m_a * 1000)))
        bundle = build_templates(m_a, n_templates, rng)
        label = format_mass_for_filename(m_a)
        dest = out_dir / f"templates_{label}.npz"
        if bundle is None:
            print(f"  m_a={m_a:.3f}: no open visible channel; skip")
            continue
        np.savez_compressed(dest, **bundle)
        brs = model.alp_branchings(m_a, model.INV_F_REF)
        print(f"  m_a={m_a:.3f}: ctau_ref={float(bundle['ctau_m_u2eq1']):.3e} m  "
              f"BR[mumu={brs['mumu']:.2f} had={brs['hadronic']:.2f} "
              f"tautau={brs['tautau']:.2f}] -> {dest.name}", flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Legacy BC10 two-track proxy templates")
    ap.add_argument(
        "--legacy-proxy",
        action="store_true",
        help="confirm generation of pre-2501 two-track proxy templates",
    )
    ap.add_argument("--mass", type=float, nargs="+", default=None)
    ap.add_argument("--n-templates", type=int, default=20_000)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--out", type=Path, default=None,
                    help="template directory (default: $GRENDEL_BC10_TEMPLATES_DIR)")
    args = ap.parse_args(argv)
    if not args.legacy_proxy:
        ap.error(
            "publication templates require: "
            "python -m grendel.models.alp.templates_pythia"
        )
    masses = args.mass if args.mass else ALP_MASS_GRID
    out_dir = args.out or ModelPaths.resolve("bc10").templates
    print(f"BC10 decay templates -> {out_dir}")
    generate(masses, args.n_templates, out_dir, args.seed)
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
