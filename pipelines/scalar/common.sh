# Shared settings of the BC4/BC5 (dark scalar) pipeline. Sourced by every step; needs `source env/activate.sh` first.
[ -n "${REPRO:-}" ] && [ -n "${GRENDEL_EXHAD_PYTHON:-}" ] || { echo "source env/activate.sh first" >&2; exit 2; }
set -euo pipefail

S="$REPRO/pipelines/scalar"                # this folder
W="$REPRO/work/scalar"                     # large data
L="$REPRO/logs/scalar"                     # logs
SRC="$REPRO/src/grendel-sensitivity"       # analysis code (8802f50; grendel/ numerically = c621279, see README)
# reference material of the validations (pipelines/common/ref.sh): REF = $GRENDEL_REF, the outputs of the original
# runs (optional: each comparison that needs it is skipped without it), the published tables and the paper figures
source "$REPRO/pipelines/common/ref.sh"
CAMP="$REF/grendel/campaigns"
EC="$CAMP/scalar/scalar_eventcalc_20261019"            # the v2.5 scans (lxplus, one job per mass)
PUB="$GRENDEL_PUBLISHED/grendel_v2_5"                  # the published tables (in this repository)
SLOT="python $REPRO/pipelines/common/slot_run.py"
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 MPLCONFIGDIR="$W/.mpl"
mkdir -p "$W" "$L" "$MPLCONFIGDIR"

# ---- the published recipe -------------------------------------------------------------------------------------
# BC4 four-vector pool (curves/data/grendel_v2/README.md: the v1 pool, scratch/bc4_nsig10_20260820):
BC4_N_POOL=400000; BC4_TILT=5; PROD_SEED=42
# produced in Aug 2026 with `particle` 0.26.2 (PDG 2024 b-hadron lifetimes); the pinned 1.0.1 gives weights
# 0.1-0.4 % lower. Only the BC4 production step runs with this overlay (00_setup.sh installs it).
PARTICLE_OVERLAY="$W/overlay/particle-0.26.2"
# BC5 four-vectors (asbuilt_20261014/chain.sh, code 461b7e5): scan_bc5 defaults (n_pool 200000, n_higgs 200000,
# seed 42, no tilt, BR(h->SS) 0.01, quartic B on), produced in ONE RNG stream over the interleaved order below.
BC5_N_POOL=200000; BC5_N_HIGGS=200000
# Higgs (pT, y) histogram: MG5 LO 200000 events per mode, Pythia 8.317 shower seed 1 (asbuilt_20261014/bc5/madgraph)
MG5_NEVENTS=200000; SHOWER_SEED=1
# exHad scalar-1809 templates (v2_baseline_20260929/templates/scalar): 20000 per mass, base seed 1234,
# sub-requests of 2000, 3 retries, chunk 512 (the templates_exhad defaults)
N_TEMPLATES=20000; TPL_SEED=1234
# scans (scalar_eventcalc_20261019/lxplus/job_sc.sh): one worker per mass, thresholds 3 10; BC5 100 decay samples
# (BC4: the default, 100)

# BC4: the 86 masses of production.MASS_GRID (= grendel_v1/bc4), in grid order (the production order)
BC4_MASSES="0.14 0.16 0.18 0.2 0.22 0.24 0.26 0.28 0.3 0.32 0.34 0.36 0.38 0.4 0.42 0.44 0.46 0.48 0.5 0.525 0.55 0.575 0.6 0.625 0.65 0.675 0.7 0.725 0.75 0.775 0.8 0.825 0.85 0.875 0.9 0.925 0.95 0.975 1 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4 4.1 4.2 4.3 4.4 4.5 4.6 4.7"
# BC5: the 112 masses in the interleaved order of asbuilt_20261014/chain.sh (the production order, which fixes the
# RNG stream) and of scalar_eventcalc_20261019/lxplus/points_bc5.txt
BC5_MASSES="0.14 0.18 0.22 0.26 0.3 0.34 0.38 0.42 0.46 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.2 2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4 4.2 4.4 4.6 4.8 5.5 7 9 11 14 16 20 25 30 40 50 60 0.16 0.2 0.24 0.28 0.32 0.36 0.4 0.44 0.48 0.525 0.575 0.625 0.675 0.725 0.775 0.825 0.875 0.925 0.975 1.05 1.15 1.25 1.35 1.45 1.55 1.65 1.75 1.85 1.95 2.1 2.3 2.5 2.7 2.9 3.1 3.3 3.5 3.7 3.9 4.1 4.3 4.5 4.7 5 6 8 10 12 15 18 22.5 27.5 35 45 55 62"

# reference data (read-only, optional: $GRENDEL_REF), used only for validation; each check is skipped without it
REF_BC5_VECTORS="$CAMP/asbuilt_20261014/bc5/llp_4vectors"           # 112 x 181 MB, readable
REF_BC4_MANIFEST="$CAMP/v2_3_20261016/lxplus/inputs_manifest.txt"   # sizes of the BC4 pool (the files are HFS-
                                                                    # compressed on the Mac drive: empty here)
REF_TEMPLATES="$CAMP/v2_baseline_20260929/templates/scalar"         # the exHad bundles of v2.5 (macOS arm64)
REF_MG5="$REF/grendel-sensitivity/grendel_work/bc5/madgraph"        # the 2026-09-25 LHE files
REF_SHOWER="$CAMP/asbuilt_20261014/bc5/madgraph"                    # their 2026-10-14 Pythia 8.317 shower

# label of a mass as in the file names (mS_0p140.csv, templates_0p140.npz)
label() { python -c "import sys; print(f'{float(sys.argv[1]):.3f}'.replace('.', 'p'))" "$1"; }

# have_ref that counts the comparisons skipped for want of reference data, for the step's OK line (step_ok)
NSKIP=0
ref_ok() { have_ref "$@" && return 0; NSKIP=$((NSKIP + 1)); return 1; }
step_ok() {
  if [ "$NSKIP" -gt 0 ]; then echo "$1 OK ($NSKIP comparison(s) skipped: reference data not available)"; else echo "$1 OK"; fi
}

# one instance of a step at a time (run_all.sh and a manual run of the same step queue up instead of
# duplicating work); the second one then finds the work done
step_lock() { exec 7>"$W/.lock.$(basename "$0" .sh)"; flock 7; }
