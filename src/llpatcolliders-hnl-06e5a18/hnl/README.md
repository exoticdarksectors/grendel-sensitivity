# `hnl/` - Heavy Neutral Lepton production and sensitivity

This package runs the full GRENDEL HNL chain for proton-proton collisions at
14 TeV:

1. produce HNL four-vectors from mesons, kaons, induced taus, prompt taus,
   and electroweak W/Z processes;
2. combine the production channels for each flavor and mass;
3. ray-cast the HNL trajectories through the GRENDEL detector geometry;
4. for the events that cross the fiducial volume, sample flavor-dependent
   FairShip decays, reconstruct the visible tracks with the shared GRENDEL
   reconstruction, and apply the PR #13 selection to get the acceptance;
5. scan the active-flavor mixing and produce sensitivity tables and plots.

The supported flavor hypotheses are `Ue`, `Umu`, and `Utau`. The committed
mass grid contains 116 points from 0.2 to 10 GeV.

## Layout

```text
hnl/
|-- README.md
|-- environment.yml
|-- config_mass_grid.py
|-- run_all.py                 production orchestrator
|-- run_analysis.py            sensitivity entry point
|-- run_full_all.sh            complete three-flavor workflow
|-- analysis/                  FairShip decay templates, GRENDEL acceptance,
|                              exclusion band, plots
|-- data/
|   `-- production/fonll/central/   committed central FONLL grids (consumed)
|-- production/
|   |-- combine_channels.py
|   |-- decay_engine/          meson, baryon, kaon, and tau decays
|   |-- fonll/                 heavy-meson spectrum parser and sampler
|   `-- madgraph/
|       |-- run_tau_production.py
|       |-- run_wz_production.py
|       `-- run_wz_sharded.py
|-- tests/
|-- vendored/
|   |-- HNLCalc/               production BRs, differential rates, lifetimes
|   |-- fairship/              FairShip HNL decay modules (rest-frame sampling)
|   `-- SM_HeavyN_CKM_AllMasses_LO/   MadGraph HeavyN model
`-- tmp/                       generated runs, caches, decay templates; git-ignored
```

The GRENDEL detector geometry and 4-hit reconstruction are not duplicated
here: the analysis imports `grendel_geometry` and `reco_common` directly from
the sibling `../higgs/` package, the single source on
exoticdarksectors/llpatcolliders `main` (PR #13).

## Physics Chain

Production includes:

- charged and semileptonic decays of `K`, `D`, `Ds`, `B`, `Bs`, and `Bc`;
- b-baryons through `Lambda_b -> Lambda_c l N` (the "Bbaryon" channel);
- induced taus from the leptonic `Ds/D+/B+/Bc -> tau nu` modes and the
  semitauonic `B/Bs -> D(*)/Ds(*) tau nu` and `Lambda_b -> Lambda_c tau nu`
  modes (polarization treatment under "Production Decay Kinematics" below);
- prompt taus from `W -> tau nu` and `gamma*/Z -> tau tau`;
- `W/Z -> ell N` through MadGraph.

The channel labels combined per mass are `Bmeson`, `Dmeson`, `Bc`,
`Bbaryon`, `tau`, `induced_tau`, `Kmeson`, and `WZ`
(`production/combine_channels.py::CHANNELS`).

The meson spectra come from the committed 14 TeV FONLL tables. HNLCalc
provides production branching ratios, differential three-body rates, and
tau-decay rates. The HNL lifetime and visible final states used downstream come
from the FairShip decay modules instead (see "Signal Decay and Acceptance").

For a production row generated at unit mixing, the downstream yield is

```text
N_signal = L_int * sum_i [
    weight_i * U_alpha^2 * P_decay_i(U_alpha^2) * acceptance_i
]
```

The analysis implements this scan at `3000 fb^-1` and uses
`N_signal >= 3` as its default threshold. Here `acceptance_i` is the GRENDEL
reconstruction + selection efficiency described next, and the lifetime
reweighting plus decay-in-volume probability `P_decay_i(U^2)` are applied in
`analysis/decay_reco_acceptance.py::scan_u2`.

### Signal Decay and Acceptance

The HNL's own decay is flavor dependent -- the visible final states differ for
`Ue`, `Umu`, and `Utau` -- so a fixed `e+e-` final state is wrong for the muon
and tau scenarios. Decays are sampled with the vendored FairShip HNL modules
(`vendored/fairship/`, driven through `ROOT.TPythia8`), which also supply the
flavor-aware lifetime `ctau(U^2 = 1)`. This is a separate template-generation
stage (`analysis/generate_decay_templates.py`) that caches rest-frame decays to
`tmp/decay_templates/<flavor>/templates_<mass>.npz`. It is the only stage that
needs PyROOT + Pythia8, and runs under the Homebrew-ROOT venv, not the conda
`hnl` env.

For each four-vector that crosses the fiducial volume, the analysis samples a
decay vertex along the flight path, boosts a FairShip template to the lab,
keeps the two highest-momentum charged stable daughters (best-two-track),
ray-casts them to the tracker walls, and runs the shared
`../higgs/reco_common` bounded 4-hit reconstruction followed by the PR #13
selection (gate / pointing / collinearity / timing). The surviving fraction is
`acceptance_i`. This is the same reconstruction and selection the cosmic-decay
background uses, so signal and background share one definition.

### Production Decay Kinematics

These are the parent-to-HNL production decays (the HNL's own decay is covered
above). Two-body decays use exact rest-frame kinematics followed by a Lorentz
boost.
The polarized two-body tau modes use a fixed unit longitudinal analyzing
power whose sign follows the tau origin: `+1` (HNL forward) for W-origin
taus, `-1` (HNL backward) for the helicity-suppressed heavy-meson leptonic
sources (`production/decay_engine/tau_decay.py`).

Meson three-body decays are sampled with
`decay_3body_weighted_dq2dE`, using the HNLCalc differential
`dBR/(dq2 dE)` expression. The b-baryon mode `Lambda_b -> Lambda_c l N` uses
`decay_3body_weighted_dq2dm122` with the HNLCalc `dBR/(dq2 dm12^2)`
expression. Leptonic tau three-body decays use
`decay_3body_weighted_dE` with the HNLCalc `dBR/dE` expression. These are not
flat phase-space samplers. The remaining approximation is the treatment of
angular and spin correlations that are not contained in those one- or
two-dimensional differential rates.

### Stored Weights

All production CSVs are headerless:

```text
weight, E, px, py, pz
```

`weight` is in pb and factors out the active mixing, so it corresponds to
`U_alpha^2 = 1`. Four-momenta are in GeV in the CMS lab frame.

Direct meson weights use

```text
2 * sigma_FONLL * fragmentation_fraction * BR(parent -> N + X) / N_sampled
```

where the factor of two converts the FONLL quark-plus-antiquark convention to
the total rate. `Bc` uses its separately configured cross section.

Induced-tau weights additionally include the parent-to-tau branching ratios
(`parent -> tau nu` for the two-body modes, `B -> D(*) tau nu` for the
semitauonic ones). The b-baryon weight carries the bottom-fragmentation closure
remainder (`FRAG_LAMBDA_B`) and `BR(Lambda_b -> Lambda_c l N)`. Prompt-tau and
W/Z rows use their MadGraph event weights and the electroweak K-factor, which
is now keyed per process (`K_FACTOR_EW_BY_PROCESS` in
`production/constants.py`: the W value for `W/Z -> ell N`, and the W vs
Drell-Yan values per tau origin) rather than a single flat constant in the
drivers. The kaon channel uses the approximate inclusive kaon flux documented
in `production/constants.py`.

## Setup

Create and activate the committed environment:

```bash
conda env create -f hnl/environment.yml
conda activate hnl
```

The direct runtime dependencies include `numpy`, `pandas`, `scipy`, `sympy`,
`mpmath`, `particle`, `matplotlib`, `tqdm`, `trimesh`, `rtree`, `numba`, and
`cycler`. `pytest` is included for tests. The full runner checks these imports
before starting expensive production.

MadGraph5_aMC@NLO v3.6.6 is required for prompt-tau and W/Z production. Its
executable is resolved in this order:

1. `$HNL_MG5_EXE`;
2. `hnl/vendored/MG5_aMC_v3_6_6/bin/mg5_aMC`;
3. the shared project-level install;
4. the legacy sibling-checkout install.

LHAPDF (binary + library) is supplied by the conda environment. Provision the
PDF data set once into the env with `lhapdf install NNPDF40_nlo_as_01180`, or
point `$HNL_LHAPDF_DATA` at an existing LHAPDF data directory. Resolution is
conda-only (no sibling-checkout fallback); a missing PDF backend fails with a
precise error before MadGraph launches.

The signal decay-template stage (`analysis/generate_decay_templates.py`)
additionally needs PyROOT + Pythia8 (`ROOT.TPythia8`), which the conda `hnl`
env does not provide. Run it under a Python that has both; on this development
machine that is the venv at
`/Volumes/sandbox/projects/aaaPHYSICSaaa/.venvs/fairship` (Homebrew ROOT with
`--with-pythia`, plus `numpy`/`scipy`). Everything else -- production,
geometry, and the exclusion scan -- runs entirely in the conda env.

## Full Run

From the repository root:

```bash
conda activate hnl
cd hnl
./run_full_all.sh
```

The script performs the complete three-flavor chain. It first runs non-W/Z
production, then runs W/Z production in independent mass shards, combines all
channels, and runs the analysis.

It does **not** generate the FairShip decay templates (that stage needs ROOT +
Pythia8, which the conda env lacks). Generate them once beforehand -- see
"Decay templates" under Manual Operation -- or the analysis will skip every
mass point with a "no decay templates" note.

By default it uses:

- 6 non-W/Z production workers;
- 12 cores for prompt-tau MadGraph generation;
- 8 simultaneous W/Z shards with one MG5 core each;
- 3 analysis workers.

Analysis workers load large combined CSVs and geometry arrays. On an 18 GB
development machine, 12 analysis workers caused severe memory pressure and a
userspace watchdog restart. Three workers completed reliably. This limit is
about memory, not CPU utilization.

The runner creates a unique tag such as `full_20260608_231500`, preventing a
new invocation from silently overwriting an older run. To resume a run, reuse
its tag:

```bash
HNL_RUN_TAG=full_20260608_231500 ./run_full_all.sh
```

W/Z production passes `--skip-existing`, so non-empty completed mass points
are retained. Other useful overrides are:

```bash
HNL_PYTHON=/path/to/python \
HNL_PRODUCTION_WORKERS=6 \
HNL_PROMPT_TAU_CORES=12 \
HNL_WZ_JOBS=8 \
HNL_WZ_CORES_PER_JOB=1 \
HNL_ANALYSIS_WORKERS=3 \
./run_full_all.sh
```

The target W/Z CPU use is `HNL_WZ_JOBS * HNL_WZ_CORES_PER_JOB`. Avoid setting
both values high enough to oversubscribe the machine.

Generated files are written below:

```text
tmp/runs/<tag>/llp_4vectors/
tmp/runs/<tag>/analysis/
tmp/cache/
tmp/wz_shards/
```

Top-level stage logs are written as `tmp/<tag>_*.log`.

## Manual Operation

### Production without W/Z

```bash
python -u run_all.py \
  --flavor Ue Umu Utau \
  --no-wz \
  --skip-combine \
  --workers 6 \
  --prompt-tau-nb-core 12
```

`run_all.py` parallelizes channel/flavor jobs. Its built-in W/Z path processes
the masses serially inside each flavor job; the full runner deliberately uses
the sharded W/Z driver instead.

### Sharded W/Z Production

```bash
python -u production/madgraph/run_wz_sharded.py \
  --flavor Umu \
  --jobs 8 \
  --nb-core 1 \
  --skip-existing
```

Run this once per flavor. Every shard has an isolated MadGraph work directory,
which avoids process-tree collisions.

### Combine

```bash
python -m production.combine_channels --flavor Ue Umu Utau
```

### Decay templates

Generate the FairShip rest-frame decay templates once per flavor (this is the
only stage that needs PyROOT + Pythia8, so use the ROOT venv, not conda):

```bash
/Volumes/sandbox/projects/aaaPHYSICSaaa/.venvs/fairship/bin/python \
  analysis/generate_decay_templates.py \
  --flavor Ue Umu Utau --n-templates 20000 --out tmp/decay_templates
```

The analysis reads `tmp/decay_templates/<flavor>/templates_<mass>.npz`; a mass
point with no template file is skipped with a note. `--skip-existing` reuses
already-generated files.

### Analysis

```bash
python -u run_analysis.py --flavor Ue Umu Utau --workers 3
```

The analysis writes:

```text
tmp/runs/<tag>/analysis/hnl_sensitivity.csv
tmp/runs/<tag>/analysis/hnl_exclusion.pdf
tmp/runs/<tag>/analysis/hnl_exclusion.png
tmp/runs/<tag>/analysis/scan_status.json
tmp/runs/<tag>/analysis/run_metadata.json
```

If an exclusion boundary is not crossed inside the configured `U^2` scan, the
CSV stores a `NaN` boundary and sets `u2_min_open` or `u2_max_open`. The plot
fills to the corresponding axis edge, omits a false closing line, and marks
the open direction.

Use `--plot-only` to regenerate the plot from an existing sensitivity CSV and
`--force-geometry` to rebuild cached ray intersections.

## Paths and Files

Path defaults are centralized in `production/paths.py`. Supported overrides
include:

```text
HNL_RUN_TAG
HNL_RUN_DIR
HNL_TMP_DIR
HNL_MG5_WORK_DIR
HNL_TAU_POOL_CSV
HNL_LLP_VECTORS_DIR
HNL_ANALYSIS_DIR
```

Mass filenames use three decimal places with `p` replacing the decimal point:
`1.025 GeV -> mN_1p025.csv`. Use
`config_mass_grid.parse_mass_from_filename` rather than parsing labels by
hand.

Zero-byte channel files are intentional threshold sentinels. They mean that
the channel is closed at that mass. Non-empty files near a threshold can have
fewer rows than the requested pool because only a subset of parent species is
kinematically allowed.

## Tests

From the repository root:

```bash
conda activate hnl
python -P -m pytest hnl/tests/ -q
```

The tests cover FONLL parsing and sampling, mass labels, channel combination,
two- and three-body kinematics, tau and kaon production, W/Z CSV conversion,
open contour handling, the FairShip->GRENDEL acceptance core (boost, best-two
tracks, reconstruction/selection, and the interior lifetime peak of the `U^2`
scan), and a meson-production smoke path. `-P` prevents any
stale ignored MG5 parser cache named `py.py` from shadowing pytest's
compatibility module.

## Physics Inputs and Limitations

- The heavy-meson backend is the committed central NNPDF4.0 NLO FONLL grid in
  `data/production/fonll/central/`, generated in the external workspace
  `/Volumes/sandbox/projects/aaaPHYSICSaaa/NNPDF40/fonll-local`.
- The FONLL tables stop at `pT = 50 GeV`. The central curve uses one central
  scale/PDF choice; scale/PDF/mass variation grids are produced in that
  external FONLL workspace and propagated into a band via
  `run_variation_band.py` (`--grid-dir`, default `tmp/fonll/output`), but are
  not folded into the central run (see `REMAINING_WORK.md`).
- Charm and bottom species share one heavy-flavor shape per table; species
  fractions are applied in the event weights. `Bc` reuses the bottom shape at
  the Bc mass. Per-species and dedicated-Bc shapes are planned
  (see `REMAINING_WORK.md`).
- b-baryons are now included through `Lambda_b -> Lambda_c l N`, which carries
  the closure remainder of the LHCb bottom-fragmentation ratios;
  `Xi_b`/`Omega_b` are not separated and are approximated as Lambda_b-like.
  The Lambda_b pT-y shape reuses the bottom FONLL grid at the Lambda_b mass.
- Bottom and charm fragmentation inputs are measured in restricted,
  different acceptances and are extrapolated as constants over the FONLL
  phase space; their uncertainties are not propagated
  (see `REMAINING_WORK.md`).
- The kaon flux is a parametrized soft-QCD approximation and dominates the
  normalization uncertainty at the lightest masses. A measured/Pythia kaon
  spectrum is planned (see `REMAINING_WORK.md`).
- Meson, baryon, and tau three-body energy/`q^2` distributions are
  HNLCalc-weighted, but complete multidimensional matrix-element spin
  correlations are not modeled (see `REMAINING_WORK.md`).
- The vendored HNLCalc form factors are not a pinned current lattice set.
  Induced-tau modes use externally normalized branching fractions, but direct
  HNL channels retain HNLCalc's absolute normalization
  (see `REMAINING_WORK.md`).
- The electroweak K-factor is a per-process table whose entries all currently
  hold the inclusive `1.3` constant; differential NLO/LO values are an optional
  upgrade (see `REMAINING_WORK.md`).
- The signal acceptance is the GRENDEL 4-hit reconstruction + PR #13 selection
  efficiency (best-two-track) applied to the FairShip decays, the same
  selection the cosmic-decay background uses. The exclusion itself is still a
  zero-background `N_signal >= 3` threshold: the cosmic background is computed
  by the shared `../higgs/` machinery but is not yet folded into the limit, and
  no systematic uncertainties are included (see `REMAINING_WORK.md`).

The FONLL and HNLCalc source details, versions, and local modifications are
recorded in `vendored/PROVENANCE.md`. Physics, detector, numerical, and
reproducibility work still needed is tracked in `REMAINING_WORK.md`.
