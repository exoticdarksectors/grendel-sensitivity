"""BC5 four-vectors exactly as the scan produces them (``grendel.models.scalar.scan_bc5.produce_missing``).

    python produce_bc5.py OUT_DIR N_POOL N_HIGGS SEED HIGGS_HIST MASS [MASS ...]

Three RNG streams (mixing, quartic B, Higgs) shared by every mass in the given order, so a file is reproduced
only by a run over the same mass list. HIGGS_HIST is the Higgs (pT, y) histogram (``tracked``: the package's
``data/higgs_pt_y_14TeV.csv.gz``, the published input). Each CSV is written to a temporary name and renamed
when complete; a rerun regenerates the whole stream and writes only the missing files.
"""
import functools
import sys
from pathlib import Path

from grendel.io.vectors import format_mass_for_filename
from grendel.models.scalar import higgs_pool, model
from grendel.models.scalar import production_bc5 as production

out = Path(sys.argv[1])
n_pool, n_higgs, seed = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
hist = sys.argv[5]
masses = [float(m) for m in sys.argv[6:]]
out.mkdir(parents=True, exist_ok=True)
dest = {m: out / f"mS_{format_mass_for_filename(m)}.csv" for m in masses}
if all(p.exists() for p in dest.values()):
    print(f"all {len(masses)} four-vector files present")
    sys.exit(0)
if hist != "tracked":
    # the production reads the tracked histogram; point its sampler at another one (no source change)
    production.sample_higgs_4vectors = functools.partial(higgs_pool.sample_higgs_4vectors, path=Path(hist))
print(f"Higgs histogram: {higgs_pool.POOL_PATH if hist == 'tracked' else hist}")

# == scan_bc5.produce_missing(masses, out, n_pool, n_higgs, seed, br_hss=0.01, quartic_b=True, force=True)
br_hss = model.BR_HSS_BC5
rngs = production.make_streams(seed)
pool, higgs, sigma_bottom = production.shared_pools(rngs, n_pool, n_higgs, br_hss=br_hss)
print(f"sigma_FONLL(bottom) = {sigma_bottom:.3e} pb; sigma(pp -> h) = {higgs['sigma_pb']:.2f} pb x "
      f"BR(h -> SS) = {br_hss:g}; {len(masses)} masses")
tmp = out / ".tmp"
tmp.mkdir(exist_ok=True)
for m in masses:
    path, counts = production.write_bc5_csv(m, tmp, rngs, pool=pool, higgs=higgs, sigma_bottom=sigma_bottom,
                                            br_hss=br_hss, quartic_b=True)
    if dest[m].exists():
        path.unlink()
        print(f"  m_S={m:.3f}: kept existing {dest[m].name}")
    else:
        path.rename(dest[m])
        print(f"  m_S={m:.3f}: mixing {counts['mixing']}, h->SS {counts['hSS']}, B->X SS {counts['BSS']} rows "
              f"-> {dest[m].name}", flush=True)
tmp.rmdir()
