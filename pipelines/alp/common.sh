# Shared settings of the BC10 (ALP) pipeline. Sourced by every step; needs `source env/activate.sh` first.
[ -n "${REPRO:-}" ] && [ -n "${GRENDEL_EXHAD_PYTHON:-}" ] || { echo "source env/activate.sh first" >&2; exit 2; }
set -euo pipefail

ALP="$REPRO/pipelines/alp"                 # this folder
W="$REPRO/work/alp"                        # large data
L="$REPRO/logs/alp"                        # logs
SRC="$REPRO/src/grendel-sensitivity"       # analysis code (8802f50; numerically = c621279 for BC10, see README)
# reference material (read only, never a production input): REF = $GRENDEL_REF, the outputs of the original runs
# (optional; without it the comparisons that need it print SKIP), the published table and the paper's PNG
source "$REPRO/pipelines/common/ref.sh"
CAMP="$REF/grendel/campaigns/bc10/bc10_v24_20261019"       # the published v2.4 run (its templates, lxplus records)
PUB="$GRENDEL_PUBLISHED/grendel_v2_4/bc10"                 # the published table (in this repository)
PAPER_PNG="$GRENDEL_PAPER_FIGURES/bc10_grendel_paper.png"  # Fig. 7 of the paper (in this repository)
SLOT="python $REPRO/pipelines/common/slot_run.py"
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 MPLCONFIGDIR="$W/.mpl"
mkdir -p "$W" "$L" "$MPLCONFIGDIR"

# The published recipe (curves/data/grendel_v2_4/README.md, bc10_v24_20261019).
N_POOL=1200000; TILT=5; PROD_SEED=42                       # production (v2.2 recipe)
TIP_RUN_MASSES="3.4 3.42 3.44 3.46 3.48 3.5 3.52 3.54 3.56 3.58"   # lxplus/job_prod.sh, one RNG stream
TIP_MASSES="3.42 3.44 3.46 3.48 3.52 3.54 3.56 3.58"       # the 8 masses taken from that run
N_TEMPLATES=20000; TPL_SEED=1234; SUB_REQUEST=250; MAX_RETRIES=10
POINTS="$ALP/points_all.txt"                               # the 153 masses of lxplus/points_all.txt
