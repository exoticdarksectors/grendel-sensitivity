#!/usr/bin/env bash
# Step 00 -- work folders and the `particle` 0.26.2 overlay for the BC4 production.
#
# The published BC4 four-vector pool (scratch/bc4_nsig10_20260820, Aug 2026) was produced with `particle` 0.26.2
# (PDG 2024); the pinned 1.0.1 of env/venv carries b-hadron lifetimes 0.06-0.4 % shorter (grendel_v2/README.md).
# 01_bc4_vectors.sh puts this overlay first on PYTHONPATH for the production only; every other step uses env/venv.
# Idempotent. Seconds.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
mkdir -p "$W"/{bc4,bc5,templates_exhad,overlay} "$L"

if ! PYTHONPATH="$PARTICLE_OVERLAY" python -c "import particle, sys; sys.exit(particle.__version__ != '0.26.2')" 2>/dev/null; then
  rm -rf "$PARTICLE_OVERLAY"
  # pure Python; its dependencies (attrs, hepunits) are those of the venv's particle 1.0.1
  python -m pip install --quiet --no-deps --target "$PARTICLE_OVERLAY" "particle==0.26.2"
fi
PYTHONPATH="$PARTICLE_OVERLAY" python - <<'EOF'
import particle
from particle import Particle
print("overlay particle", particle.__version__, "tau(B+) =", Particle.from_pdgid(521).lifetime, "ns")
EOF
python -c "import particle; from particle import Particle; print('venv    particle', particle.__version__, 'tau(B+) =', Particle.from_pdgid(521).lifetime, 'ns')"
echo "00_setup OK"
