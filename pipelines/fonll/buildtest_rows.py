"""Build-flag test: a few raw quark-grid rapidity rows with one FONLL build (compare with buildtest_cmp.py).

    source env/activate.sh && source pipelines/fonll/env.sh
    FONLL_OPT="" FONLL_BUILD_DIR=work/fonll/variants/O0 bash pipelines/fonll/01_build_fonll.sh       # a variant
    python pipelines/common/slot_run.py -n 1 -- python pipelines/fonll/buildtest_rows.py work/fonll/fonll-1.3.3 bottom work/fonll/variant_test/main 0,50
    python pipelines/fonll/buildtest_cmp.py work/fonll/variant_test/main/nl_b_.out $FONLL_ORIG_WS/run/nlo_bottom/nl_b_.out 0,50

Arguments: FONLL tree, quark, output dir, comma-separated indices into generate.Y_VALUES[, LHAPDF id].
Results (2026-11-06) in README.md, Build notes.
"""
import os, subprocess, sys, time
from pathlib import Path
from grendel.production.fonll_grids import generate as gen

build, quark, outdir = Path(sys.argv[1]).resolve(), sys.argv[2], Path(sys.argv[3]).resolve()
iys = [int(x) for x in sys.argv[4].split(",")]
pdfid = int(sys.argv[5]) if len(sys.argv) > 5 else 331700
outdir.mkdir(parents=True, exist_ok=True)
ys = [gen.Y_VALUES[i] for i in iys]
mass = gen.QUARKS[quark]["mass"]
prefix = f"nl_{quark[0]}_"
card = gen.grid_input(prefix, pdfid, mass, y_values=ys)
for f in outdir.glob(prefix + "*"):
    f.unlink()
t0 = time.time()
with open(outdir / "fonllgrid.log", "w") as log:
    subprocess.run([str(build / "Linux" / "fonllgridlha")], cwd=outdir, input=card, text=True, stdout=log,
                   stderr=subprocess.STDOUT, check=True)
print(f"{build.name} {quark} rows {iys}: {time.time() - t0:.1f} s")
