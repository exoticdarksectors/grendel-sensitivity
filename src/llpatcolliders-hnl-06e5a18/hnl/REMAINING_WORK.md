# Remaining work

This is the authoritative inventory of work still needed before the generated
GRENDEL HNL exclusion curves can be presented as publication-grade physics
results. It is derived from the code paths that produce the curves, not from
comments or earlier planning notes.

The inventory covers four different kinds of work:

- **P0 -- result definition:** missing detector, background, statistical, or
  physics-model inputs that can change the meaning of the exclusion.
- **P1 -- uncertainty/completeness:** missing channels and uncertainty
  propagation needed for defensible central curves and bands.
- **P2 -- validation/reproducibility:** numerical convergence, provenance, and
  workflow checks needed to make the result repeatable and auditable.

An item is not complete merely because a central constant exists in the code.
It is complete when its source, uncertainty, propagation, and regression test
are present. The induced-tau central branching fractions in
`production/constants.py` have been manually reviewed against PDG, lattice,
and LHCb inputs, but that review is not yet reproducible from a machine-readable
source/version table. Their provenance, uncertainties, and correlations remain
outstanding below.

A companion diagnostic estimate of the curve impact is stored in
`audits/curve_impact_20260610/README.md`, with detailed results at every
configured mass point and in 0.3 GeV bins. Those ranges are explicitly
classified as calculations, literature-informed ranges, proxies, or unknowable
scenarios; they are not a statistical uncertainty band.

## Curve-visibility triage (production x analysis, by PBC scenario)

The items below are ordered by physics priority, not by their effect on the
*current* curves. As of 2026-06-16 every PBC scenario -- 100 (Ue), 010 (Umu),
001 (Utau) -- has complete production (116 combined CSVs each) and complete
decay templates (116 each); the live selection in
`analysis/decay_reco_acceptance.py` matches the upstream `higgs/` GRENDEL
reconstruction cut-for-cut, and the exclusion threshold is the standard
zero-background 95% CL value (`N >= 3`). No remaining item moves any of the three
curves visibly: each is an NLO refinement (sub-leading uncertainty, provenance,
convergence) or a model gap that changes the *meaning* of the result
(backgrounds, detector response) rather than the plotted central line.

| layer \ scenario | 100 (Ue) | 010 (Umu) | 001 (Utau) |
|---|---|---|---|
| **production** | mesons + W/Z EW; kaon flux parametric (NLO) | mesons + W/Z EW; kaon flux parametric (NLO) | mesons + tau-parent + W/Z EW; tau Kallen & W-tau polarization already fixed |
| **analysis** | cuts/geometry/threshold single-sourced with `higgs/` -- aligned | aligned | aligned |

Robustness fixes that *prevent* a future silent curve corruption (not a current
defect) were applied on 2026-06-16: the geometry cache now invalidates when its
source CSV is newer (mtime guard, item 19); `run_sensitivity.py` records every
skipped `(flavor, mass)` with its reason in `run_metadata.json` and exits
non-zero when nothing is processed; and `run_full_all.sh` generates decay
templates before the analysis so a clean checkout cannot emit an empty plot. The
dead, stale `SEP_MIN/SEP_MAX/P_CUT` block was removed from
`analysis/constants.py` (the live cuts were always sourced from
`decay_reco_acceptance.py`).

## P0: define the detector-level result

### 1. Add detector response to the decay/reconstruction acceptance

**Current code:** the old analytic two-body acceptance (`analysis/sensitivity.py`)
is gone. `analysis/decay_reco_acceptance.py` now runs a per-event Monte Carlo:
for each HNL four-vector that crosses the fiducial volume it samples decay
vertices along the flight path, draws flavor-dependent FairShip rest-frame decay
templates, boosts them to the lab, keeps the two highest-momentum charged stable
daughters (best-two-track), and runs the shared `../higgs/reco_common` bounded
4-hit reconstruction plus the PR #13 selection (gate/pointing/collinearity/
timing). The visible branching fraction is implicit (the fraction of templates
with two reconstructable charged tracks), not a single inclusive `BR_vis`
factor. What is still idealized: detector response is geometric wall hits with
Gaussian position/time smearing only -- no material interactions, tracking or
vertexing inefficiency, trigger/readout, pileup, dead regions, or occupancy --
and the FairShip templates are Pythia-level final states with no detector
simulation.

**Required work:**

- implement detector response beyond geometric hits + smearing: tracking,
  vertexing, particle thresholds, material interactions, trigger/readout, and
  event-selection efficiencies;
- include pileup, timing, dead regions, and occupancy if they affect the
  proposed detector;
- assign per-channel / per-flavor efficiency uncertainties and validate the
  FairShip visible-final-state modelling against an independent generator.

**Completion test:** an efficiency map or detector simulation with versioned
inputs reproduces benchmark samples, and the sensitivity code consumes
per-channel efficiencies with uncertainty variations.

### 2. Validate the geometry and detector configuration

**Current code:** the shared `../higgs/grendel_geometry.py` (the PR #13 single
source) builds a closed tunnel mesh from one survey polyline, a fixed 22 m
vertical position, and a fixed 24 cm wall inset, and classifies tracker vs
scintillator surfaces. Rays originate at `(0, 0, 0)`. The mesh contains no
supports, services, inactive regions, material, alignment uncertainty, or
configurable detector layout.

**Required work:**

- confirm the survey coordinates, CMS IP transform, tunnel cross-section, and
  proposed active volume with the detector/design owners;
- version the geometry and expose layout parameters rather than relying on
  import-time constants;
- model inactive and inaccessible regions and evaluate alignment/survey
  variations;
- document the coordinate convention and independently compare ray
  intersections against a reference geometry implementation.

**Completion test:** signed-off geometry inputs, reference intersection tests,
and a geometry-variation envelope propagated to the curves.

### 3. Supply backgrounds and a statistical model

**Current code:** `analysis/constants.py` defines exclusion as
`N_signal >= 3`, a zero-background 95% CL approximation. Signal now passes the
PR #13 selection (gate/pointing/collinearity/timing), and a cosmic
decay-in-flight background is modelled by the shared `../higgs/` machinery, but
that background is not folded into the limit: there is still no background count
in the statistical model, no control region, nuisance parameter, systematic
uncertainty, or coverage calculation.

**Required work:**

- estimate beam-, collision-, cosmic-, neutrino-, and detector-induced
  backgrounds after the final selection;
- define control samples and background uncertainties;
- include luminosity, reconstruction, trigger, geometry, signal-model, and
  finite-MC nuisance parameters;
- choose and implement the statistical prescription (for example CLs or a
  documented Poisson construction) and validate its expected coverage;
- define how open upper/lower sensitivity boundaries are reported.

**Completion test:** the final likelihood/configuration is versioned, has toy
or asymptotic validation, and produces expected limits with uncertainty bands.

### 4. Fix the physics-hypothesis definition

**Current code:** production and lifetime tables are evaluated for the three
single-flavor patterns `Ue`, `Umu`, and `Utau`. The exclusion labels only
`|U|^2`; the Majorana/Dirac convention, charge-conjugate counting, and the
relationship between production mixing and total lifetime are not stated in
the result metadata.

**Required work:**

- state and validate the Majorana/Dirac convention used by HNLCalc and the
  MadGraph UFO;
- audit all factors of two and charge-conjugate channels under that convention;
- support arbitrary mixing ratios if the intended result is not restricted to
  the three single-flavor hypotheses;
- regenerate lifetime, visible branching fraction, production, and decay
  inputs consistently for every advertised hypothesis.

**Completion test:** each output curve carries machine-readable model metadata
and agrees with independent benchmark rates for the same convention.

## P1: production model and its uncertainties

### 5. Produce FONLL scale, PDF, heavy-quark-mass, and grid-coverage variations

**Current code:** the committed bottom and charm grids contain one central
scale choice and the central `NNPDF40_nlo_as_01180` member. They stop at
`pT = 50 GeV` and `|y| = 3`. No perturbative, PDF, heavy-quark-mass, alpha_s,
or phase-space truncation uncertainty reaches the exclusion curves.

**Required work:**

- standard seven-point `(mu_R, mu_F)` variations;
- all NNPDF4.0 NLO replicas, combined with the documented Monte Carlo
  prescription;
- bottom- and charm-quark mass variations consistent with the FONLL setup;
- optional alpha_s variations if they are not covered by the selected PDF set;
- extended `pT` and rapidity grids, or a quantitative bound on the omitted
  contribution after GRENDEL acceptance;
- retain both normalization and shape variations rather than reducing every
  variation to a single total-cross-section factor.

Use the same three-column rectangular table format as the committed grids:

```text
# pT_GeV  rapidity  d2sigma_dpTdy_pb_per_GeV
```

Add a variation manifest containing the FONLL revision, PDF ID/member,
`m_b`/`m_c`, scales, beam energy, grid bounds, and checksums.

**Completion test:** central and varied grids can be regenerated from the
manifest, and the full production plus analysis chain yields an uncertainty
band at every mass point.

### 6. Replace shared heavy-flavor shapes and static fragmentation fractions

**Current code:** `B+`, `B0`, and `Bs` share one bottom shape; `D0`, `D+`, and
`Ds` share one charm shape. Species masses are rebuilt after sampling.
`FRAG_B` and `FRAG_C` are constant fractions measured in restricted and
different acceptances, extrapolated over the full grid without covariance.
Fragmentation-function parameters and feeddown uncertainties are absent.

**Required work:**

- provide species-specific `d2sigma/dpTdy` grids or differential
  fragmentation fractions for each weakly decaying species;
- propagate the experimental fraction covariance and its `pT`/rapidity
  dependence;
- vary the fragmentation-function model and parameters used in the FONLL
  calculation;
- document strong-decay feeddown conventions so no parent is double counted.

**Completion test:** species fractions close consistently at each kinematic
point, reproduce their source measurements in the quoted acceptance, and
produce correlated curve variations.

### 7. Add omitted charmed-baryon HNL sources

**Current code:** `OMITTED_FRAG_C` records approximately 36.7% of the charm
fragmentation allocation as `Lambda_c+`, `Xi_c0`, `Xi_c+`, and `J/psi`, but
none contributes to direct HNL production or to the induced-tau pool.

**Required work:**

- identify and implement the relevant `Lambda_c` and `Xi_c` leptonic or
  semileptonic HNL production channels with current form factors;
- provide their production spectra and fragmentation fractions with
  uncertainties;
- evaluate whether charmonium HNL modes are relevant over the 0.2--10 GeV
  mass grid and either implement or quantitatively dismiss them;
- audit any additional weakly decaying charm species omitted from the current
  closure.

**Completion test:** the charm fragmentation accounting is explicit and every
omitted component has either a generated channel or a documented negligible
bound on the final curves.

### 8. Replace the inclusive bottom-baryon closure approximation

**Current code:** the entire bottom-fragmentation remainder is treated as
`Lambda_b`, with the bottom FONLL shape rebuilt at the `Lambda_b` mass and only
`Lambda_b -> Lambda_c l N` generated. `Xi_b` and `Omega_b` are not separated.

**Required work:**

- provide differential production fractions and spectra for `Lambda_b`,
  `Xi_b`, and `Omega_b`;
- implement their relevant HNL and tau-producing decays with current form
  factors;
- propagate the large baryon-fraction kinematic dependence and covariance.

**Completion test:** bottom species close without an unidentified remainder,
and alternative baryon models are included in the production uncertainty.

### 9. Replace the borrowed Bc model

**Current code:** `Bc` uses the bottom FONLL shape at the `Bc` mass and an
unvalidated inclusive normalization `SIGMA_BC_PB = 0.9e6 pb`. That same model
also normalizes `Bc -> tau nu`.

**Required work:**

- generate a dedicated 14 TeV `Bc+ + Bc-` spectrum using BcVegPy, NRQCD, or a
  documented equivalent;
- provide normalization, scale, PDF, heavy-quark-mass, and NRQCD/model
  uncertainties;
- keep the charge convention explicit;
- propagate the correlated normalization to direct-Bc and induced-tau
  channels.

**Completion test:** a versioned grid and uncertainty set replace both the
shape reuse and `SIGMA_BC_PB`, with independent cross-section benchmarks.

### 10. Replace the kaon production and transport model

**Current code:** `generate_kaon_csvs.py` samples only charged kaons using a
Tsallis `pT` model and a Gaussian rapidity with
`SIGMA_KAON_PB = 3.0e11 pb`. Neutral `K_S/K_L -> pi l N` modes available in
HNLCalc are omitted. Every kaon is decayed immediately, only the HNL momentum
is stored, and analysis rays every HNL from IP5. Parent lifetime, magnetic
bending, material survival, interaction losses, and the displaced HNL
production vertex are therefore absent. This model can dominate the
lowest-mass result.

**Required work:**

- produce charged- and neutral-kaon spectra from measured data and/or tuned
  Pythia 8 over the full kinematic region that can reach GRENDEL;
- state whether the input is a per-event yield or a cross section and provide
  the corresponding inelastic normalization;
- vary generator tune, soft-QCD model, multiplicity, charge ratio, and
  extrapolation outside measured coverage;
- propagate each kaon through the magnetic field and detector material,
  sample its decay or interaction point, and originate the HNL there;
- include `K_S` and `K_L` lifetimes, regeneration/interactions, and branching
  channels where relevant;
- include secondary kaon production/absorption only if it belongs in the
  signal definition, without double counting generator particles.

**Completion test:** data/generator comparisons, transport validation, and an
uncertainty envelope replace the prompt charged-kaon approximation in
kinematics, origin, and normalization.

### 11. Separate and upgrade electroweak production

**Current code:** LO MadGraph samples use one NNPDF central member and a flat
`K = 1.3`. Direct `W -> l N` and `Z -> nu N` rows are merged without a parent
tag, so all rows receive the W factor. Prompt-tau W and Drell-Yan origins are
separable but use the same central factor. Because the W and Drell-Yan factors
are currently identical (`K = 1.3`), the missing per-row W/Z split has **no
numerical effect on the present curves**; it is a logic refinement that only
matters once the factors differ.

**Required work:**

- carry the LHE mother/process identity into direct-HNL CSV rows;
- calculate process-, mass-, and preferably kinematics-dependent NLO QCD and
  electroweak corrections for W, Z, and Drell-Yan tau production;
- propagate MG scale/PDF variations and matching/shower uncertainties if a
  showered sample is used;
- check process definitions, interference, off-shell contributions, and
  overlap with heavy-flavor or top-origin W samples;
- implement or quantitatively dismiss Higgs-mediated `h -> nu N` production
  for every advertised model hypothesis.

**Completion test:** W, Z, and Drell-Yan yields can be varied independently
and reproduce reference inclusive and differential cross sections.

### 12. Close the heavy-hadron production-channel inventory

**Current code:** direct meson production is a hand-written list of exclusive
two- and three-body modes. It is not checked against the larger HNLCalc channel
inventory; for example, HNLCalc contains `D+ -> omega l N` while the production
list omits it. Higher-resonance and nonresonant semileptonic modes are not
accounted for. Weak cascades such as `B -> D/Ds X` followed by charm-hadron HNL
production are also absent.

**Required work:**

- generate a channel-closure table comparing implemented modes with all
  kinematically allowed modes in the selected rate model;
- add material missing exclusive, higher-resonance, and nonresonant channels,
  or bound their effect;
- include weak heavy-hadron cascades with decay vertices and avoid overlap
  with prompt heavy-flavor production;
- perform the same closure audit for tau-producing parent decays beyond the
  selected 11 sources.

**Completion test:** every allowed parent fraction and decay class is
implemented or assigned a documented negligible bound, with cascade
double-counting tests.

## P1: decay, lifetime, and polarization model

### 13. Replace or pin the HNLCalc rate model

**Current code:** `vendored/HNLCalc` has no retained upstream commit or release
tag. Direct meson and baryon channels use its embedded form factors and
absolute branching ratios. Several parameterizations are old or weakly
sourced. Manual `m_N -> 0` spot checks found substantial differences between
selected HNLCalc integrals and external tau-mode branching fractions, but the
comparison script and benchmark table are not yet committed. Two concrete
defects have already been found and fixed in the vendored copy (a misplaced
Källén-factor parenthesis in `tau -> P N` and an `|U|^4` mixing factor in the
baryon channel; see `vendored/PROVENANCE.md`), both also present upstream --
which strengthens the case for a full channel-by-channel validation.

**Required work:**

- pin an upstream HNLCalc revision or replace it with a maintained local rate
  library;
- update decay constants, CKM inputs, hadron masses, and semileptonic form
  factors to a coherent current set;
- include form-factor covariance and normalization uncertainties for every
  material direct channel, including baryons and Bc;
- validate integrated rates and differential shapes against independent
  calculations and SM-limit branching fractions for all three flavors;
- validate and assign uncertainty to the abrupt exclusive-hadron versus
  inclusive-parton width matching around 1 GeV, including continuity and
  quark-hadron-duality assumptions;
- make the direct-HNL and induced-tau channels use compatible rate inputs.

**Completion test:** channel-by-channel benchmark tables and regression tests
cover integrated rates and selected differential distributions across mass.

### 14. Propagate parent-to-tau branching-fraction uncertainties

**Current code:** the 11 induced-tau sources use manually reviewed central
values for `Ds`, `D+`, `B+`, `B`, `Bs`, `Bc`, and `Lambda_b` decays. The review
is not encoded as a reproducible source/version table, and no
experimental/theory uncertainty or covariance is propagated. The `Bs` and
`Bc` entries are theory-derived, and the `Lambda_b` measurement contains an
external-normalization uncertainty.

**Required work:**

- store source, date/version, uncertainty components, and correlations beside
  every central value;
- propagate HFLAV/PDG experimental averages and lattice/theory uncertainties;
- correlate common inputs such as `Vcb`, decay constants, normalization modes,
  and form factors;
- define an update procedure rather than manually changing constants.

**Completion test:** a machine-readable input table generates the constants
and produces correlated induced-tau normalization variations.

### 15. Validate the FairShip lifetime and visible-fraction model

**Current code:** the HNL lifetime `ctau(U^2 = 1)` and the visible final states
now come from the FairShip decay templates (`analysis/generate_decay_templates.py`):
the lifetime is `HNLbranchings.computeNLifetime`, and the visible fraction is
implicit in the per-flavor template multiplicity. The previous
`data/ctau/*.dat` lifetime / `BR_vis` tables and the `BR_vis = 1` fallback are
gone. Still unrecorded: the FairShip module revision, the exact decay-mode and
visible definition for GRENDEL, the Majorana/Dirac convention, finite-template
statistics, and any uncertainty on the lifetime or visible fraction. Templates
are generated on a fixed mass grid and matched per mass label (no interpolation
between masses).

**Required work:**

- pin the FairShip module revision and write a manifest (mass grid,
  `n_templates`, seed, couplings, decay-selection config);
- define exactly which modes count as visible for GRENDEL and confirm the
  template-derived visible fraction matches that definition;
- propagate lifetime and visible-fraction uncertainties and correlate them with
  production;
- quantify finite-template statistics and confirm the mass-grid spacing needs
  no interpolation (or add it);
- cross-check the FairShip total width / lifetime and visible fraction against
  an independent calculation.

**Completion test:** a pinned, manifested template generation, strict loading,
threshold-aware mass coverage, and regression benchmarks against an independent
lifetime/visible-fraction calculation.

### 16. Complete tau decay modes and implement realistic spin correlations

**Current code:** semileptonic energy and invariant-mass distributions are
HNLCalc-weighted, but remaining angles are sampled isotropically. Two-body
tau-to-HNL decays use fixed unit-analyzing-power longitudinal asymmetries
with the origin-dependent sign (`+1` for W-origin taus, `-1` for the
helicity-suppressed heavy-meson leptonic sources; see
`production/decay_engine/tau_decay.py`). The magnitude is still an
approximation: the true analyzing power is mass-dependent and below one for
the vector modes, the asymmetry is applied about the lab tau direction
(no Wigner rotation, which matters most for `Ds -> tau nu` where the tau is
nearly at rest in the parent frame), and tau charge is not retained in the
prompt-tau CSV. Semitauonic, baryonic, and Drell-Yan tau sources are treated
as unpolarized, although their polarization depends on charge, phase space,
and production kinematics. Tau-to-HNL decay generation includes only `pi`,
`K`, `rho`, `K*`, and leptonic three-body modes; multi-hadron spectral modes
such as `a1/3pi` are absent.

**Required work:**

- use full matrix-element event generation or helicity amplitudes for the
  relevant parent-to-tau-to-HNL and direct-HNL chains;
- retain tau spin density matrices through production and decay;
- include Drell-Yan and W tau polarization with their kinematic dependence;
- add all material multi-hadron tau decay modes using validated spectral
  functions or a spin-aware generator such as TAUOLA;
- demonstrate branching-fraction closure as a function of HNL mass;
- validate angular and energy distributions against an independent generator.

**Completion test:** spin-correlated benchmark distributions and a quantified
curve shift relative to the current approximation.

## P1/P2: numerical and workflow uncertainty

### 17. Quantify finite Monte Carlo and integration errors

**Current code:** most sampled pools default to 100,000 rows, HNLCalc
three-body integrals commonly use 500 samples, and rejection-sampling ceilings
are numerical approximations (the weighted samplers now detect a proposal
exceeding the ceiling and restart unbiased with a raised ceiling, so an
under-estimated grid scan can no longer silently truncate the density, but
grid resolution still sets the proposal efficiency). The curve contains no
finite-sample uncertainty.

**Required work:**

- run independent seeds and increasing pool sizes by channel and mass;
- propagate weighted-event variance through geometry and the U-squared scan;
- validate rejection envelopes and retry limits near thresholds;
- converge the HNLCalc integrations and differential samplers;
- allocate events adaptively to rare sources and low-acceptance regions.

**Completion test:** per-mass MC/integration uncertainty is below a declared
target or appears in the final band.

### 18. Converge the acceptance, serialization, and exclusion scans

**Current code:** the accepted decay probability is a Monte Carlo estimate over
`DECAY_SAMPLES = 100` decay vertices sampled uniformly along the in-volume path,
with the reconstruction + selection re-evaluated per vertex
(`analysis/decay_reco_acceptance.py::build_event_mc`/`scan_u2`). The mixing scan
has 200 points from `1e-12` to `1e-1`, followed by local interpolation. The mass
grid has a minimum spacing of 15 MeV and becomes much coarser at high mass;
plotted lines simply connect calculated points.

**Required work:**

- increase or replace the eight-significant-digit CSV serialization and
  quantify its effect on reconstructed mass, boost, geometry, and acceptance,
  especially for ultra-boosted low-mass HNLs;
- demonstrate the decay-vertex Monte Carlo is converged (vary `DECAY_SAMPLES`
  and seeds) against an exact quadrature of the decay-in-volume integral;
- replace or validate the fixed mixing scan with bracketed root finding for
  both exclusion boundaries;
- demonstrate stability under mass-grid refinement, especially at production
  thresholds and where the dominant channel changes;
- do not claim sub-grid precision: the current grid cannot support estimates
  every 0.3 MeV.

**Completion test:** curve boundaries are stable under all three refinements
within a declared numerical tolerance.

### 19. Make caches and channel completeness self-validating

**Current code:** geometry caches are named only by flavor and mass and do not
yet embed a checksum of the input vectors, geometry version, origin, or cuts.
Since 2026-06-16 a cache is invalidated when its source 4-vector CSV is newer
than the cached NPZ (mtime guard), so regenerating a combined CSV forces a fresh
ray-cast; a full content/geometry checksum is still outstanding (use
`--force-geometry` after geometry-code edits). `run_sensitivity.py` now records
every skipped `(flavor, mass)` with its reason in `run_metadata.json` and exits
non-zero when no point is processed, so dropped masses are no longer silent.
Channel combination is strict by default everywhere since 2026-06-11:
`combine_channels.py` fails on missing or malformed channels while accepting
zero-byte closed-channel sentinels, partial runs must opt in via
`--allow-missing`, and `run_full_all.sh` therefore fails on an incomplete run.

**Required work:**

- key caches by input checksum and geometry/analysis configuration;
- write a run manifest with code revision, environment, seeds, input
  checksums, channel row counts, and summed weights;
- verify that no stale channel or cache from another run tag is consumed.

**Completion test:** deliberately modifying an input invalidates the cache,
and incomplete runs fail before producing a combined curve.

### 20. Pin the software and external tools

**Current code:** the Python environment, MadGraph, LHAPDF, PDF data,
compiler/runtime, and some vendored sources are not captured in one immutable
run description. The HNLCalc snapshot is specifically unpinned.

**Required work:**

- commit an environment lock and record exact external executable versions;
- checksum vendored tables, UFO files, and generated process cards;
- preserve MG run cards, param cards, logs, and random seeds for final samples;
- add a one-command clean-room reproduction of representative production and
  analysis points.

**Completion test:** a fresh checkout can reproduce benchmark weights and
curve points within the stated stochastic tolerance.

### 21. Expand regression coverage for the new production paths

**Current code:** `tests/test_production_channels.py` (2026-06-11) covers the
11-source induced-tau table (source count, polarization tags, kinematic
thresholds), the induced-tau and `Bbaryon` drivers through their production
interfaces (on-shell HNLs, weight normalization, closed-mass sentinels, the
`dq2dm122` sampler), and `tests/test_combine_channels.py` distinguishes
missing channels from zero-byte sentinels under the strict default. The
stale future-dated test metadata has been removed.

**Required work:**

- tighten the 30% normalization tolerances once item 13 provides
  deterministic benchmark rates (HNLCalc's integrator is Monte Carlo at
  `nsample=500`, so the current checks pin factor-2/|U|^4-class bugs only).

### 22. Propagate all variations into final plots and tables

**Current code:** `plot_exclusion.py` draws only central connected curves. It
does not show theory, detector, background, numerical, or luminosity bands,
nor does it overlay external constraints. Run metadata does not identify all
inputs used.

**Required work:**

- define a variation registry and correlated combination procedure;
- save central, individual variation, and combined-envelope boundaries at
  every calculated mass;
- show open scan boundaries honestly and avoid interpolation across missing or
  insensitive masses;
- include input/version metadata and, if required for the intended figure,
  current external constraints with source/version tracking;
- publish machine-readable result tables in addition to PDF/PNG figures.

**Completion test:** every visible band can be traced to stored variation
outputs and regenerated from the run manifest.

## Work that can be done in the NNPDF40 workspace

The active FONLL grid-generation workspace is external to this repository:
`/Volumes/sandbox/projects/aaaPHYSICSaaa/NNPDF40/fonll-local`. Do not maintain
a second active copy under this HNL tree; use the external workspace for item 5
and part of item 6. Status as of 2026-06-11:

1. **Done.** `scripts/generate_meson_grids.py --campaign` accepts
   `(mu_R, mu_F)` (verified against FONLL's `read ffact,fren` order), PDF
   members, and `m_b`/`m_c` variations, with per-grid tags, headers, and a
   SHA-256 manifest. `validate_variation_gate.py` is the hard-abort
   central-reproduction gate.
2. **Pilot done.** A members 0-10 + 7-point scale + mass campaign (38 grids)
   ran on 2026-06-11; the regenerated centrals reproduce the committed
   grids. The full 100-replica campaign (218 grids, multi-day) is still
   outstanding.
3. **Done with caveat.** `combine_variations.py` writes scale/PDF/mass/
   combined envelopes plus a manifest. The PDF statistics use replica
   members >= 1 only (member 0 excluded; an earlier revision included it and
   biased the band low). Envelope grids are pointwise constructions marked
   `envelope_band` and are refused by `production/fonll/fonll_parser.py`;
   only coherent individual grids may be sampled.
4. **Known artifact.** The charm `(0.5, 0.5)` scale point drives `mu_F`
   below the NNPDF4.0 grid minimum at low pT (suppression x14-20 for
   pT <~ 2 GeV); see `fonll-local/PUBLICATION_BASELINE.md` for the
   prescription caveat before quoting the low-mass charm scale band.
5. **Outstanding:** the alpha_s companion grids (the `NNPDF40_nlo_as_01170`/
   `_01190` sets are installed with all 101 members and wired into
   `generate_meson_grids.py` as `--pdf nlo_as_01170` / `--pdf nlo_as_01190`;
   the four central grids still have to be generated), extended
   `pT`/rapidity coverage to bound the GRENDEL-accepted tail (item 6), and
   per-species fragmentation outputs/variations (item 8).

A full replica campaign takes days and produces large logs unless run with
`--compress-logs`; it is resumable via `--reuse-existing-grids`.

## Recommended execution order

1. Define the detector signal, geometry, backgrounds, and statistical model
   (items 1--4).
2. Pin/regenerate HNL rates, lifetime, and visible fractions (items 13--15).
3. Produce FONLL uncertainty grids and fix the largest production models:
   kaons, charm baryons, Bc, fragmentation, and channel closure (items 5--12).
4. Upgrade HNL/tau decay completeness and spin correlations (items 13--16).
5. Complete numerical convergence, strict provenance, regression coverage,
   and final uncertainty propagation (items 17--22).

Until the P0 items are resolved, curve-shift estimates for the remaining items
are diagnostics of the current simplified analysis, not uncertainties on a
fully defined experimental exclusion.
