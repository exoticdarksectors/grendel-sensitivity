#!/bin/bash
# HNL step 00: the production venv work/hnl/venv_prod with the package versions of the published HNL production
# (the Mac "hnl" conda env, logs/mac_hnl_conda_freeze.txt of campaigns/hnl/hnl_fine_condor_20261014; the lxplus
# venv_prod of the fine grid copied them). particle 0.26.2 matters: HNLCalc takes the parent lifetimes from it,
# and particle 1.0.1 (env/venv) has b-hadron lifetimes 0.06-0.4 % shorter, i.e. different production weights.
# The scan runs in env/venv (Mac scan venv: numpy 2.5.3, not available for Python 3.11; see README).
# The complete package set (direct and transitive, exact versions) is pipelines/hnl/venv_prod.lock.txt, the pip
# freeze of the venv the reproduction ran in; the direct pins of the published production are kept below.
# Resume-safe: does nothing when the venv already has exactly the locked versions.
set -euo pipefail
source "$(dirname "$0")/../../env/activate.sh"
V=$REPRO/work/hnl/venv_prod
LOCK=$REPRO/pipelines/hnl/venv_prod.lock.txt
# direct pins (the Mac "hnl" env); sympy 1.14.0 as on lxplus (it pins mpmath<1.4, the Mac env has 1.4.1)
#   numpy==2.4.6 scipy==1.17.1 pandas==3.0.3 particle==0.26.2 hepunits==2.4.5 numba==0.65.1 llvmlite==0.47.0
#   mpmath==1.4.1 matplotlib==3.10.9 cycler==0.12.1 six==1.17.0 tqdm==4.68.0 sympy==1.14.0
locked() { grep -v '^#' $LOCK; }
if [ -x $V/bin/python ] && diff <($V/bin/pip freeze --exclude grendel-sensitivity 2>/dev/null) <(locked) > /dev/null; then
  echo "venv_prod ok: $V"; exit 0
fi
"$REPRO/env/conda/bin/python3.11" -m venv --clear $V
$V/bin/pip install -q --upgrade pip
$V/bin/pip install -q --no-deps -r <(locked)
# the analysis package from src/ (no install: PYTHONPATH keeps the venv free of a second copy)
echo "$REPRO/src/grendel-sensitivity" > $V/lib/python3.11/site-packages/grendel_src.pth
$V/bin/pip freeze > $REPRO/work/hnl/venv_prod.freeze.txt
diff <($V/bin/pip freeze --exclude grendel-sensitivity) <(locked) || { echo "venv_prod differs from $LOCK" >&2; exit 1; }
$V/bin/python -c "import grendel, particle, numpy; print('venv_prod ready:', grendel.__file__, 'particle', particle.__version__, 'numpy', numpy.__version__)"
