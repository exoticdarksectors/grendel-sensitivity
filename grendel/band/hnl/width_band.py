"""Seam-derived HNL hadronic-width uncertainty band ``delta(m)``.

``delta(m)`` is the fractional uncertainty on the HNL total width
``Gamma_tot`` from quark-hadron duality, read off FairShip's *own*
``max(mesonWidth, quarkWidth)`` disagreement (the matching seam in
``NDecayWidth``):

* **floor** (``FLOOR``) below the quark threshold ``MTHR`` -- the HNL decays
  exclusively to ``pi/K/rho`` + leptons there, all measured, so the raw seam
  (where ``quarkWidth = 0``) is meaningless and is floored;
* ``|quarkWidth - mesonWidth| / Gamma_tot`` through the duality transition --
  the genuine matching ambiguity where the exclusive sum starts missing
  multi-meson states and the tree-level quark width is not yet reliable;
* **capped** at ``CAP`` above -- the inclusive quark width is trusted to its
  perturbative/duality error (~20%, Bondarenko et al. arXiv:1805.08567); the
  raw seam keeps growing only because the *exclusive* meson list is
  incomplete, which is not an uncertainty.

The raw ``|q-m|`` metric has a *structural zero* where ``quarkWidth ==
mesonWidth`` (the ~1.2-1.4 GeV crossover) -- precisely the most
duality-uncertain point -- so it is **enveloped**: through the transition
``[MTHR, m_full]`` (``m_full`` = the mass where the raw seam first reaches
``CAP``) the band is ``max(raw seam, linear ramp FLOOR->CAP)``, so the
crossover does not notch the band down through the region that should be
widest.

``CAP = 20%`` is deliberately conservative for a projection upper edge: the
width model already applies a perturbative ``QCD_correction``, so the honest
*residual* duality error is ~10%.

The same ``delta(m)`` drives the decay-model band coherently through both the
lifetime leg (``ctau = hbar/Gamma_tot``) and the decay-composition leg
(``vis_frac``); see ``decay_model_band``.

Computing the widths needs FairShip (PyROOT), so the table is precomputed once
in a ROOT-enabled environment with FairShip fetched by ``third_party/fetch.py``::

    python -m grendel.band.hnl.width_band     # writes models/hnl/data/width_band_delta.csv

and the analysis loads the CSV at scan time (no ROOT needed).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from ...io.paths import third_party_dir

TABLE_PATH = Path(__file__).resolve().parents[2] / "models" / "hnl" / "data" / "width_band_delta.csv"

# 5% exclusive floor, 20% QCD/duality cap, quark threshold ~1 GeV.
FLOOR, CAP, MTHR = 0.05, 0.20, 1.0
FLAVORS = ("Ue", "Umu", "Utau")
_COUPLINGS = {"Ue": [1.0, 0.0, 0.0], "Umu": [0.0, 1.0, 0.0], "Utau": [0.0, 0.0, 1.0]}


def delta_from_widths(mass: float, meson_w: float, quark_w: float, gtot: float) -> float:
    """Fractional Gamma_tot uncertainty from the meson/quark width seam."""
    if gtot <= 0 or mass < MTHR:
        return FLOOR
    return float(np.clip(abs(quark_w - meson_w) / gtot, FLOOR, CAP))


# ---- table generation (needs ROOT/FairShip) --------------------------------

def _seam_widths(mass: float, flavor: str):
    """(mesonWidth, quarkWidth, Gamma_tot) in GeV from FairShip HNLbranchings."""
    fairship = third_party_dir() / "FairShip" / "python"
    if str(fairship) not in sys.path:
        sys.path.insert(0, str(fairship))
    import hnl as fhnl  # requires PyROOT
    h = fhnl.HNLbranchings(float(mass), _COUPLINGS[flavor])
    lep = h.Width_3nu() + h.Width_charged_leptons()
    mes = h.Width_neutral_mesons() + h.Width_charged_mesons()
    qrk = h.Width_quarks_neutrino() + h.Width_quarks_lepton()
    return mes, qrk, lep + max(mes, qrk)


def _transition_envelope(masses, raw_seam):
    """Lift the band through the duality transition so the |q-m| structural
    zero at the meson/quark crossover does not notch it down: ramp FLOOR->CAP
    from MTHR to ``m_full`` and take ``max(raw seam, ramp)``; below MTHR stay
    at FLOOR."""
    masses = np.asarray(masses, float)
    raw = np.asarray(raw_seam, float)
    above = masses >= MTHR
    hit = np.where(above & (raw >= CAP))[0]
    m_full = float(masses[hit[0]]) if len(hit) else float(masses[above].max())
    ramp = FLOOR + (CAP - FLOOR) * (masses - MTHR) / max(m_full - MTHR, 1e-6)
    out = np.full(len(masses), FLOOR)
    out[above] = np.clip(np.maximum(raw[above], ramp[above]), FLOOR, CAP)
    return out


def compute_table(masses=None) -> pd.DataFrame:
    if masses is None:
        masses = np.round(np.concatenate([np.arange(0.20, 1.0, 0.05),
                                           np.arange(1.0, 6.01, 0.10)]), 3)
    frames = []
    for flavor in FLAVORS:
        recs = []
        for m in masses:
            mes, qrk, gtot = _seam_widths(float(m), flavor)
            raw = abs(qrk - mes) / gtot if gtot > 0 else 0.0
            recs.append({"flavor": flavor, "mass_GeV": float(m),
                         "meson_over_gtot": mes / gtot if gtot > 0 else 0.0,
                         "quark_over_gtot": qrk / gtot if gtot > 0 else 0.0,
                         "gtot_GeV": gtot,
                         "seam_raw": raw,
                         "delta_floor_only": delta_from_widths(float(m), mes, qrk, gtot)})
        df = pd.DataFrame(recs)
        df["delta"] = _transition_envelope(df["mass_GeV"], df["seam_raw"])
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


# ---- table consumption (no ROOT) -------------------------------------------

_CACHE = None


def load_delta_table() -> pd.DataFrame:
    global _CACHE
    if _CACHE is None:
        if not TABLE_PATH.exists():
            raise FileNotFoundError(
                f"{TABLE_PATH} missing; run `python -m grendel.band.hnl.width_band` "
                "in a ROOT-enabled environment to precompute it.")
        _CACHE = pd.read_csv(TABLE_PATH)
    return _CACHE


def delta(flavor: str, mass: float) -> float:
    """delta(m) for one (flavor, mass), interpolated from the precomputed table."""
    t = load_delta_table()
    sub = t[t["flavor"] == flavor].sort_values("mass_GeV")
    return float(np.interp(mass, sub["mass_GeV"].to_numpy(), sub["delta"].to_numpy()))


def delta_and_hadfrac(flavor: str, mass: float):
    """(delta, had_frac) where had_frac = Gamma_had/Gamma_tot = max(meson, quark)
    fraction. The composition leg needs had_frac to convert a delta on Gamma_tot
    into the hadronic-channel reweight (Gamma_had +/- delta*Gamma_tot)."""
    t = load_delta_table()
    sub = t[t["flavor"] == flavor].sort_values("mass_GeV")
    m = sub["mass_GeV"].to_numpy()
    d = float(np.interp(mass, m, sub["delta"].to_numpy()))
    hf = float(np.interp(mass, m, np.maximum(
        sub["meson_over_gtot"].to_numpy(), sub["quark_over_gtot"].to_numpy())))
    return d, hf


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Recompute the seam-derived width-uncertainty table delta(m) "
                                             "from FairShip's HNL widths (needs PyROOT and third_party/FairShip).")
    ap.add_argument("--out", type=Path, default=TABLE_PATH, help=f"output CSV (default: {TABLE_PATH})")
    args = ap.parse_args(argv)
    df = compute_table()
    df.to_csv(args.out, index=False)
    n = len(df) // len(FLAVORS)
    print(f"wrote {args.out} ({len(df)} rows, {n} masses x {len(FLAVORS)} flavors)")
    for fl in FLAVORS:
        s = df[df["flavor"] == fl]
        print(f"  {fl}: delta range [{s['delta'].min():.3f}, {s['delta'].max():.3f}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
