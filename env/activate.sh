# Source this file (bash) before any step of the reproduction:
#     source env/activate.sh
# It puts the user-space toolchain (env/conda), the analysis venv (env/venv) and the
# Wolfram Engine (env/wolfram, if installed) on PATH and points the grendel package at
# the folders of this checkout. Nothing outside this directory is used.

REPRO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export REPRO
# Reference material of the validation steps (read only; no production step uses it): the published tables and the
# paper's figure files are in this repository; the outputs of the original runs (GRENDEL_REF) are optional, see
# reference/README.md. Without it the comparisons that need it are skipped. Default: ../extra_space, if it exists
# (the machine the reproduction was made on); GRENDEL_REF=none runs without it.
export GRENDEL_PUBLISHED="$REPRO/src/grendel-sensitivity/curves/data"
export GRENDEL_PAPER_FIGURES="$REPRO/reference/paper_figures"
if [ -z "${GRENDEL_REF:-}" ]; then
  if [ -d "$(dirname "$REPRO")/extra_space/grendel" ]; then GRENDEL_REF="$(dirname "$REPRO")/extra_space"; else GRENDEL_REF=none; fi
fi
export GRENDEL_REF

export PATH="$REPRO/env/venv/bin:$REPRO/env/conda/bin:$REPRO/env/bin:$PATH"
[ -d "$REPRO/env/wolfram/bin" ] && export PATH="$REPRO/env/wolfram/bin:$PATH"

# compilers of the conda toolchain (FONLL, MadGraph, Pythia)
export CC="$REPRO/env/conda/bin/x86_64-conda-linux-gnu-gcc"
export CXX="$REPRO/env/conda/bin/x86_64-conda-linux-gnu-g++"
export FC="$REPRO/env/conda/bin/x86_64-conda-linux-gnu-gfortran"
# env/conda/lib is deliberately NOT put on LD_LIBRARY_PATH (it breaks the system curl/git);
# conda-built binaries carry their rpath. For a program that cannot find libLHAPDF at run time:
#     LD_LIBRARY_PATH="$REPRO/env/conda/lib" <program>
export GRENDEL_CONDA_LIB="$REPRO/env/conda/lib"

# external software (third_party/fetch.py checks out here)
export GRENDEL_THIRD_PARTY_DIR="$REPRO/third_party"
export GRENDEL_LHAPDF_CONFIG="$REPRO/env/conda/bin/lhapdf-config"
export LHAPDF_DATA_PATH="$REPRO/work/lhapdf"
export GRENDEL_LHAPDF_DATA="$LHAPDF_DATA_PATH"
export GRENDEL_MG5_EXE="$REPRO/third_party/MG5_aMC_v3_6_6/bin/mg5_aMC"
export GRENDEL_PYTHIA8_CONFIG="$REPRO/third_party/pythia8317/bin/pythia8-config"
export PYTHIA8DATA="$REPRO/third_party/pythia8317/share/Pythia8/xmldoc"
export GRENDEL_EXHAD_PYTHON="$REPRO/third_party/exHad/.venv/bin/python"

# large data and results stay under work/ (git-ignored)
export GRENDEL_WORK_DIR="$REPRO/work/grendel_work"
export GRENDEL_RESULTS_DIR="$REPRO/work/results"

# a fixed date stamp in the PDFs (matplotlib's /CreationDate), so that rerendering a figure gives the same bytes:
# 2026-11-08, the day the last pipeline of the reproduction finished
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1794096000}"

# one thread per process: parallelism is by process (one job per mass)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
