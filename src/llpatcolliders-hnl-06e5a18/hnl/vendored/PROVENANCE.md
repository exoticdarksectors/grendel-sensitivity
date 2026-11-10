# Vendored sources

## FONLL meson tables

The default production backend uses two local FONLL+LHAPDF
double-differential heavy-meson production cross sections, `dsigma/dpT/dy`,
generated for the HL-LHC pp 14 TeV setup.

- Generator release: [`leoredi/fonll-nnpdf40` v0.1.0](https://github.com/leoredi/fonll-nnpdf40/releases/tag/v0.1.0),
  commit `34cd8b8`. The vendored `.dat` files here are byte-identical to
  the v0.1.0 release artifacts except for one header comment line
  (`# meson_mass: ...`), which is annotated in the PR1 copies to document
  the `ifrframe = 1` / `xmh` insensitivity argument captured below.
  Regeneration: clone the upstream tag, apply the two patches against
  FONLL v1.3.3, run
  `python scripts/generate_meson_grids.py --pdf nlo --quark {bottom,charm}`.
- FONLL source: v1.3.3 (`c7086e49141cf6705cf7a4bc5f7d0b3a38673203`)
- LHAPDF: 6.5.6
- Process: pp at sqrt(s) = 14 TeV
- PDF set: `NNPDF40_nlo_as_01180` (LHAPDF ID 331700)
- Output: meson-level, fragmentation fraction 1
- Bottom convention: public FONLL B-hadron default, Kartvelishvili
  `alpha = 24.2`
- Charm convention: public FONLL `D0`, BCFY `r = 0.1`, with calibrated
  `D* -> D0` feeddown
- Grid: 100 pT nodes over 0..50 GeV, 100 y nodes over -3..3
- Scale: central (mu_R = mu_F = mu_0)
- Fragmentation frame: `ifrframe = 1` (y=0-frame fragmentation, public-FONLL
  default for `dsigma/dpT/dy`). Under this convention `xsecfrag` itself
  does not reference `xmh`; the only `xmh` dependence in `fragmfonll.f`
  on the `ifrframe = 1` branch is the kinematic guard at
  `szmin` (`fragmfonll.f:648`), `(pT^2 + xmh^2) * cosh(y)^2 > sh/4`.
  In the phase space of these tables (`pT <= 50 GeV`, `|y| <= 3`,
  `sh/4 = 4.9e7 GeV^2`) the guard's worst case
  `(50^2 + 5.28^2) * cosh(3)^2 ~ 2.5e5 GeV^2`
  is always satisfied, so the table is insensitive to `xmh` here —
  verified empirically by regenerating with `xmh = M_B0 / M_D0` and
  recovering byte-identical numerics. The sampler in
  `meson_sampler.py` (`sample_meson_4vectors`, and the shared-shape helper
  `meson_4vec_from_kinematics`) reconstructs the on-shell meson four-vector
  with the physical PDG mass, which is the kinematically consistent
  thing to do for this convention. If the grid is ever extended toward
  the kinematic edge of `sh/4`, the guard would become active and this
  insensitivity claim would need to be re-verified.
- Local FONLL capacity patch: `fonllgrid.f` and `fragmfonll.f` rapidity
  array limits raised to 120 nodes so the dense grid and terminating sentinel
  are accepted; no physics formulas changed
- Files:
  - `fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_charm.dat`
  - `fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_bottom.dat`

Species splitting (D0/D+/Ds for charm; B+/B0/Bs for bottom) is applied
downstream by `hnl/production/constants.py::FRAG_C` and `FRAG_B`. The
FONLL tables supply pT-y shape and normalization before physical species
fractions.

## MadGraph HeavyN UFO model (W/Z -> ell N production)

`SM_HeavyN_CKM_AllMasses_LO/` is the FeynRules/UFO model used by the
electroweak production path (`hnl/production/madgraph/run_wz_production.py`).

- Model: HeavyN, CKM / AllMasses LO variant.
- Reference: Degrande, Mattelaer, Pascoli, Ruiz et al., arXiv:1602.06957
  (updated by Pascoli et al., arXiv:1812.08750).
- Provenance: copied verbatim from the upstream `llpatcolliders_FONLL`
  vendored MadGraph install
  (`vendored/MG5_aMC_v3_6_6/models/SM_HeavyN_CKM_AllMasses_LO`); only the
  build artifact `__pycache__/` was stripped.
- The driver loads it with an absolute `import model <path>` so MadGraph
  picks up this vendored copy regardless of which MG5 install is used.

### MadGraph binary (not vendored here)

The MadGraph5_aMC@NLO v3.6.6 install (~148 MB) is intentionally **not**
committed. `production/madgraph/_mg5_common.py::_resolve_mg5_exe` resolves
the `mg5_aMC` executable at runtime in this order:

1. `$HNL_MG5_EXE` (explicit path), else
2. `hnl/vendored/MG5_aMC_v3_6_6/bin/mg5_aMC` (drop-in vendoring if you want
   a self-contained checkout), else
3. the projects-root shared `vendored/MG5_aMC_v3_6_6/bin/mg5_aMC` (two
   levels above `hnl/`, shared across sibling `llpatcolliders_*` checkouts), else
4. the sibling `llpatcolliders_FONLL/vendored/MG5_aMC_v3_6_6/bin/mg5_aMC`.

To make the checkout self-contained, copy the upstream install into
`hnl/vendored/MG5_aMC_v3_6_6/` (consider git-lfs or a `.gitignore` entry for
the binary tree given its size).

## HNLCalc

Pure-Python computation of HNL production and decay branching ratios.

- Upstream: https://github.com/laroccod/HNLCalc (Feng, Hewitt, Kling,
  La Rocco; arXiv:2405.07330).
- Vendored from our fork https://github.com/leoredi/HNLCalc, branch
  `physics-fixes` at commit `e292cef9` — i.e. upstream `laroccod/HNLCalc`
  `main` @ `07f84728` plus the four local-fix commits listed under "Local
  modifications" below. `HNLCalc.py` here is byte-identical to that fork
  commit (normalized to LF line endings; the prior snapshot was checked out
  CRLF).
- Vendored files: `HNLCalc.py`, `alph_str.csv`, `README.md`
- Stripped from the snapshot: `.git/`, notebook examples, embedded PNG plots,
  build artifacts (`__pycache__/`).
- Form-factor caveat: most embedded form-factor parameterizations lack
  channel-level citations, so this snapshot is not independently traceable as
  a current lattice/HFLAV form-factor set. The SM-limit normalization audit is
  documented in `REMAINING_WORK.md`.

Used in `hnl/production/decay_engine/` to evaluate:
- `get_2body_br` and `get_3body_dbr_*` for meson -> HNL channels
- `get_3body_dbr_baryon` for the b-baryon channel
  (`Lambda_b -> Lambda_c l N`)
- `get_2body_br_tau` and `get_3body_dbr_tau` for tau -> HNL channels
- `integrate_3body_br` for numerical phase-space integration; the b-baryon
  channel calls it with `integration='dq2dm122'`

Production calls HNLCalc with unit coupling (`U^2 = 1`) so the active mixing
factors cleanly from each stored row. `analysis/run_sensitivity.py` performs
the explicit `U^2` scan; the HNL lifetime and visible final states it reweights
with come from the FairShip decay templates (below), not from HNLCalc. HNLCalc
itself is used only for production BRs, three-body differential rates, tau-decay
rates, and phase-space integration (the calls listed above).

### Local modifications

These are the four commits on `leoredi/HNLCalc@physics-fixes` (pinned above).

- Commented out the unused `from skhep.math.vectors import LorentzVector,
  Vector3D` (HNLCalc.py line ~15). Those symbols are never referenced in the
  module; the import only forced a hard dependency on scikit-hep's
  `skhep.math` (dropped from modern scikit-hep), so the file would not import
  standalone. No numerical effect.

- Removed a duplicated `Ds+ -> K0` form-factor block in `HNLCalc.py`.
  Upstream contained two independent `if` statements matching the same
  condition (`pid0 in ["431","-431"] and pid1 in ["311","-311"]`), each
  assigning `f00` a different value: the first set `f00 = 0.747` and a
  later block set `f00 = 0.72`. Because both are plain `if`s, the second
  always overrode the first, making the `0.747` block dead code. We
  deleted the dead `0.747` block and kept the effective value
  (`f00 = 0.72`), so the numerical behavior is unchanged.

  `Ds+ -> K0` is a c -> d transition (with an anti-strange spectator),
  distinct from the `D -> K` c -> s transition immediately above it in
  `HNLCalc.py`. The retained normalization is supported by Melikhov and
  Stech (hep-ph/0001113, Table XVII), which gives
  `F_+^{Ds -> K}(0) = F_0^{Ds -> K}(0) = 0.72` in a relativistic
  dispersion constituent-quark model. For comparison, HPQCD gives
  `f_+^{D -> K}(0) = 0.747(19)` for `D -> K` (1008.4562), while an
  independent lattice calculation gives `f_+^{Ds -> K}(0) = 0.68(4)(3)`
  (0903.1664), compatible with `0.72` within its quoted uncertainties.
  The deleted `0.747` block appears to have copied the `D -> K` value
  into the `Ds -> K` channel.

- Fixed a misplaced parenthesis in `get_2body_br_tau` (pseudoscalar
  branch, `tau -> P N` for `P = pi, K`). Upstream builds the two-body
  phase-space (Källén) factor as `sqrt(1 - A*(1 - B))` instead of
  `sqrt((1 - A)*(1 - B))`, with `A = ((M1 - m_N)/m_tau)^2` and
  `B = ((M1 + m_N)/m_tau)^2`; the vector branch three lines above has the
  correct parenthesization. The upstream form overestimates
  `BR(tau -> pi N)` by x1.07 at `m_N = 0.5 GeV`, x1.38 at `1.0 GeV`,
  x3.9 at `1.5 GeV`, and x8.5 at `1.6 GeV` (`tau -> K N` similarly:
  x1.9 at `1.0 GeV`). The same line is present in the public upstream
  repository (github.com/laroccod/HNLCalc) as of 2026-06; this fix is a
  deliberate local deviation.

- Fixed a double-squared mixing factor in `get_3body_dbr_baryon`.
  Upstream sets `Ulx = vcoupling^2` and then multiplies the squared
  amplitude by `Ulx^2`, giving `|U|^4` instead of `|U|^2`. Changed `Ulx`
  to the unsquared coupling so the `Ulx^2` in the amplitude yields
  `|U|^2`. Numerically inert for this package's single-flavor instances
  (`vcoupling` is 0 or 1), but wrong for any mixed-coupling use.

## FairShip HNL decay modules (`fairship/`)

`fairship/` holds the FairShip (github.com/ShipSoft/FairShip) Heavy-Neutral-
Lepton physics modules, used to sample **flavor-dependent** HNL rest-frame
decays (the visible final states differ for `Ue`/`Umu`/`Utau`) and to supply
the flavor-aware lifetime `ctau(U^2 = 1)`. They are unmodified copies imported
via the project's `hnl_alaship` package. Per-file roles, the decay-selection
config (`analysis/fairship_decay_selection.conf`), and the PyROOT + Pythia8
template-generation flow are documented in `fairship/PROVENANCE.md`.

## Prompt-tau pool (`tau_pool.csv`)

`tmp/cache/tau_pool.csv` is the default LHE-extracted tau four-vector pool
used by the prompt-tau driver (`production/madgraph/run_tau_production.py`,
Stage 1). The older `vendored/tau_pool.csv` remains readable as a
compatibility fallback, but new Stage-1 generation writes to `hnl/tmp`.
The pool is a generated data product: a flavor- and m_N-independent input
that downstream loops consume many times.

### What it is

- Single MG5 v3.6.6 run of the proc card
  `production/madgraph/cards/proc_card_tau_production.dat`:
    - `pp -> W+ -> tau+ nu_tau`
    - `pp -> W- -> tau- nu_tau_bar`
    - `pp -> tau+ tau-`   (direct Drell-Yan; covers gamma*/Z exchange and
                           gamma*/Z interference)
- LHE parsed by `production/madgraph/lhe_to_csv.py::write_tau_csv`, which
  emits a 6-column headerless CSV: `w, E, px, py, pz, origin`.
- `origin` is the mother PDG resolved via MOTHUP1 in the LHE event block:
  `+/-24` for W-mediated, `23` for explicit Z propagator, anything else
  (gluon, light quarks) for direct-DY events.
- The per-row `w` already encodes `sigma_LO / N_pool_events`. For
  `tau+ tau-`, each tau row keeps the full event weight because either tau is
  an independent potential HNL parent. The pool is pre-multiplied by the
  process-keyed EW K-factor per row (`K_FACTOR_EW_BY_PROCESS["W"]` for
  W-origin taus, `["DY"]` for Drell-Yan taus; both currently `1.3`, see
  `production/constants.py`) so the summed weight approximates `sigma_NLO`.

### MG5 / PDF settings

- MG5_aMC v3.6.6, SM model (`import model sm`); no HeavyN at this stage.
- `define p = g u c d s b u~ c~ d~ s~ b~` (5-flavor proton, see proc card).
- PDF: `pdlabel = lhapdf`, `lhaid = 331700`
  (`NNPDF40_nlo_as_01180`, same set the FONLL meson tables use).
- LHAPDF is discovered at runtime from the active environment or the fallback
  locations implemented in `production/madgraph/_mg5_common.py`.

### Regenerate

```
cd hnl
rm -f tmp/cache/tau_pool.csv
python -m production.madgraph.run_tau_production \
    --nevents 100000 --nb-core 4
```

(`--skip-mg5` is a `store_true` flag — omit it to regenerate, pass it to
reuse the cached pool. Or simply delete the file and re-run
`python run_all.py`; the pipeline will rebuild Stage 1 automatically.)
