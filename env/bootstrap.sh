#!/usr/bin/env bash
# Build everything the pipelines need, in user space, on a generic x86_64 Linux box
# (no root, no system packages beyond bash, curl, git, tar, bzip2 and a C library).
#
#     bash env/bootstrap.sh            # idempotent: skips what is already there
#     bash env/bootstrap.sh --unlocked # resolve the newest versions instead of the lock files
#     bash env/bootstrap.sh --minimal  # steps 1-3 only: enough for figures/redraw.sh (minutes, ~2 GB)
#
# 1. micromamba 2.10.0 (static binary) -> env/bin/micromamba
# 2. conda toolchain env/conda (python 3.11, gcc/g++/gfortran, cmake, LHAPDF 6.5): the exact builds of
#    env/conda-linux-64.lock (--unlocked: env/environment.yml)
# 3. analysis venv env/venv: the exact versions of env/requirements-venv.lock.txt, then
#    src/grendel-sensitivity (editable). The lock holds the extras test,production,raycast,digitize,
#    HNLimits 1.2.0 without its dependencies (they would pull a matplotlib that garbles the log
#    ticks) and six (MadGraph).
# 4. the external software at the pins of src/grendel-sensitivity/third_party/PINS.json
# 5. the PDF set NNPDF40_nlo_as_01180 -> work/lhapdf
# 6. Pythia 8.317 (plain build, as exHad requires), MadGraph first-run configuration
# 7. exHad in its own venv (env/requirements-exhad.lock.txt), built against that Pythia
# The Wolfram Engine (SensCalc, BC10 widths) is installed by pipelines/senscalc; it needs a
# (free) Wolfram ID to activate.
set -euo pipefail
REPRO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPRO"
mkdir -p logs work
LOCKED=1 MINIMAL=0
for a in "$@"; do
  case "$a" in
    --unlocked) LOCKED=0 ;;
    --minimal) MINIMAL=1 ;;
    *) echo "unknown option $a" >&2; exit 2 ;;
  esac
done

if [ ! -x env/bin/micromamba ]; then
  mkdir -p env/bin
  curl -fsSL https://micro.mamba.pm/api/micromamba/linux-64/2.10.0 | tar -xj -C env bin/micromamba
fi
if [ ! -x env/conda/bin/python3.11 ]; then
  if [ $LOCKED = 1 ]; then CONDA_SPEC=env/conda-linux-64.lock; else CONDA_SPEC=env/environment.yml; fi
  MAMBA_ROOT_PREFIX="$REPRO/env/mamba" env/bin/micromamba create -y -p "$REPRO/env/conda" -f "$CONDA_SPEC"
fi
if [ ! -x env/venv/bin/python ]; then
  env/conda/bin/python3.11 -m venv env/venv
fi
# shellcheck disable=SC1091
source env/activate.sh
if [ $LOCKED = 1 ]; then
  pip install -q pip==26.3
  pip install -q --no-deps -r env/requirements-venv.lock.txt
  pip install -q --no-deps -e src/grendel-sensitivity
else
  pip install -q --upgrade pip
  pip install -q -e "src/grendel-sensitivity[test,production,raycast,digitize]"
  pip install -q --no-deps HNLimits==1.2.0
  pip install -q six
fi
if [ $MINIMAL = 1 ]; then echo "bootstrap done (minimal: analysis venv only)"; exit 0; fi

python src/grendel-sensitivity/third_party/fetch.py exHad HNLCalc SM_HeavyN_CKM_AllMasses_LO fonll MG5_aMC_v3_6_6 pythia8317
python src/grendel-sensitivity/third_party/fetch.py --verify || true

mkdir -p "$LHAPDF_DATA_PATH"
for f in lhapdf.conf pdfsets.index; do
  [ -e "$LHAPDF_DATA_PATH/$f" ] || cp "env/conda/share/LHAPDF/$f" "$LHAPDF_DATA_PATH/"
done
if [ ! -f "$LHAPDF_DATA_PATH/NNPDF40_nlo_as_01180/NNPDF40_nlo_as_01180.info" ]; then
  # the set the reproduction used (DataVersion 1); the tarball is checked against its sha256 of 2026-11-10
  curl -fsSL -o work/NNPDF40_nlo_as_01180.tar.gz https://lhapdfsets.web.cern.ch/current/NNPDF40_nlo_as_01180.tar.gz
  echo "a2ef486909853a5a23ea65bb7e9b71e98425e5af12a5106d1450d0b7bb5d290d  work/NNPDF40_nlo_as_01180.tar.gz" | sha256sum -c -
  tar -xzf work/NNPDF40_nlo_as_01180.tar.gz -C "$LHAPDF_DATA_PATH" && rm work/NNPDF40_nlo_as_01180.tar.gz
fi

if [ ! -f third_party/pythia8317/lib/libpythia8.so ]; then
  (cd third_party/pythia8317 && ./configure --prefix="$PWD" && make -j"$(nproc)")
fi

if [ ! -f third_party/MG5_aMC_v3_6_6/input/mg5_configuration.txt ]; then
  printf 'set auto_update 0\nset automatic_html_opening False\nsave options\nexit\n' > work/mg5_init.txt
  (cd third_party/MG5_aMC_v3_6_6 && python bin/mg5_aMC "$REPRO/work/mg5_init.txt")
fi

if [ ! -x third_party/exHad/.venv/bin/python ]; then
  env/conda/bin/python3.11 -m venv third_party/exHad/.venv
fi
if ! third_party/exHad/.venv/bin/python -c 'import exhad' 2>/dev/null; then
  if [ $LOCKED = 1 ]; then
    third_party/exHad/.venv/bin/pip install -q --no-deps -r env/requirements-exhad.lock.txt
    EXHAD_PIP="-q --no-deps -e ."
  else
    EXHAD_PIP="-q -e ."
  fi
  # shellcheck disable=SC2086
  (cd third_party/exHad && .venv/bin/pip install $EXHAD_PIP \
     && .venv/bin/python tools/setup.py --pythia8-dir "$REPRO/third_party/pythia8317" --jobs "$(nproc)")
fi
echo "bootstrap done"
