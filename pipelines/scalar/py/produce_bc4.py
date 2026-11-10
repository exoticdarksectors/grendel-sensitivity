"""BC4 four-vectors exactly as the scan produces them (``grendel.models.scalar.scan.produce_missing``).

    python produce_bc4.py OUT_DIR N_POOL SEED TILT MASS [MASS ...]

One RNG stream: the b-hadron pool, then every mass in the given order, so a file is reproduced only by a
run over the same mass list. Each CSV is written to a temporary name and renamed when complete, so an
interrupted run leaves no truncated file; a rerun regenerates the whole stream (the production is minutes)
and rewrites only missing files: files present are kept as they are (no rewrite, mtime unchanged).
"""
import sys
from pathlib import Path

import numpy as np
import particle

from grendel.models.scalar import production
from grendel.production.fonll.fonll_parser import get_sigma_total
from grendel.production.fonll.meson_sampler import sample_meson_4vectors

out = Path(sys.argv[1])
n_pool, seed, tilt = int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
masses = [float(m) for m in sys.argv[5:]]
out.mkdir(parents=True, exist_ok=True)
print(f"particle {particle.__version__} ({particle.__file__})")
dest = {m: out / f"mS_{production._mass_label(m)}.csv" for m in masses}
if all(p.exists() for p in dest.values()):
    print(f"all {len(masses)} four-vector files present")
    sys.exit(0)

# == scan.produce_missing(masses, out, n_pool, seed, force=True, high_pt_tilt_scale=tilt), written atomically
rng = np.random.default_rng(seed)
sigma_bottom = get_sigma_total("bottom")
pool = sample_meson_4vectors(n_pool, "bottom", rng=rng, high_pt_tilt_scale=tilt, nominal_mixture_fraction=0.5)
print(f"sigma_FONLL(bottom) = {sigma_bottom:.3e} pb; {len(masses)} masses, pool {n_pool}, tilt {tilt}, seed {seed}")
tmp = out / ".tmp"
tmp.mkdir(exist_ok=True)
for m in masses:
    # the draws of every mass are made, written or not, to keep the stream
    path, n = production.write_scalar_csv(m, tmp, n_pool, rng, sigma_bottom=sigma_bottom, pool=pool)
    if dest[m].exists():
        path.unlink()
        print(f"  m_S={m:.3f}: kept existing {dest[m].name}")
    else:
        path.rename(dest[m])
        print(f"  m_S={m:.3f}: {n} rows -> {dest[m].name}", flush=True)
tmp.rmdir()
