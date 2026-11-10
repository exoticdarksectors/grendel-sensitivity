"""One rapidity chunk of a central raw quark grid, or the merge of all chunks (step 03).

    python pipelines/fonll/raw_chunk.py --quark bottom --chunk 2 --nchunks 4     # one fonllgridlha process
    python pipelines/fonll/raw_chunk.py --quark bottom --merge --nchunks 4       # -> run/nlo_bottom_central/nl_b_.out

This is grendel.production.fonll_grids.generate.build_quark_grid with --grid-workers N, split so that every
chunk is its own process and slot (the shipped driver runs the N chunks as threads of one process, and an
interrupted run restarts all of them): the same input card (gen.grid_input with the chunk's rapidities from
gen._contiguous_chunks(gen.Y_VALUES, N)), the same prefix and run directory layout, the same merge (header of
the first chunk, the bodies in rapidity order, 1 + 100 x 80 lines). The merged file is the one
`generate --reuse-existing-grids` then fragments (gen.raw_grid_complete). Chunking does not change the numbers
(byte-identical 4-chunk and 1-process raw grids on macOS, fonll-local run/nlo_bottom vs run/nlo_bottom_central).

Resume-safe: a complete chunk is kept; an incomplete one is deleted and recomputed (fonllgridlha appends).
"""
import argparse
import sys
import time
from pathlib import Path

from grendel.production.fonll_grids import generate as gen

ap = argparse.ArgumentParser()
ap.add_argument("--quark", required=True, choices=sorted(gen.QUARKS))
ap.add_argument("--nchunks", type=int, default=4)
ap.add_argument("--chunk", type=int)
ap.add_argument("--merge", action="store_true")
ap.add_argument("--check-only", action="store_true", help="exit 0 if the chunk (or the merged grid) is complete, else 1")
args = ap.parse_args()

ws = gen.Workspace.resolve()                       # GRENDEL_FONLL_DIR / GRENDEL_FONLL_GRID_DIR (pipelines/fonll/env.sh)
var = gen.central_variation(args.quark)
pdf_key = "nlo"
pdf_id = gen.PDFS[pdf_key][1] + var.pdf_member
run_dir = ws.run / f"{pdf_key}_{args.quark}_{var.tag}"   # as generate()
prefix = f"{pdf_key[:2]}_{args.quark[:1]}_"
grid_file = run_dir / f"{prefix}.out"
chunks = gen._contiguous_chunks(gen.Y_VALUES, args.nchunks)
chunk_root = run_dir / f"grid_chunks_n{args.nchunks}"


def chunk_file(i: int) -> Path:
    return chunk_root / f"{i:02d}" / f"{prefix}.out"


def chunk_complete(i: int) -> bool:
    f = chunk_file(i)
    if not f.exists():
        return False
    text = f.read_text()
    return text.endswith("\n") and len(text.splitlines()) == 1 + len(chunks[i]) * gen.FONLL_GRID_NPT


if args.merge:
    if gen.raw_grid_complete(grid_file):
        print(f"{grid_file} complete")
        sys.exit(0)
    missing = [i for i in range(len(chunks)) if not chunk_complete(i)]
    if missing:
        print(f"chunks not complete yet: {missing}")
        sys.exit(3)
    combined, header = [], None
    for i in range(len(chunks)):
        lines = chunk_file(i).read_text().splitlines()
        if header is None:
            header = lines[0]
            combined.append(header)
        elif lines[0] != header:
            raise SystemExit(f"inconsistent raw-grid header in {chunk_file(i)}")
        combined.extend(lines[1:])
    if len(combined) != 1 + len(gen.Y_VALUES) * gen.FONLL_GRID_NPT:
        raise SystemExit(f"merged raw grid has {len(combined)} lines")
    tmp = grid_file.with_suffix(".out.part")
    tmp.write_text("\n".join(combined) + "\n")
    tmp.rename(grid_file)
    print(f"merged {len(chunks)} chunks -> {grid_file}")
    sys.exit(0)

i = args.chunk
if i is None or not 0 <= i < len(chunks):
    raise SystemExit(f"--chunk must be in 0..{len(chunks) - 1}")
if gen.raw_grid_complete(grid_file) or chunk_complete(i):
    print(f"[{args.quark} chunk {i}] complete")
    sys.exit(0)
if args.check_only:
    sys.exit(1)
d = chunk_root / f"{i:02d}"
d.mkdir(parents=True, exist_ok=True)
for f in d.glob(f"{prefix}*"):
    f.unlink()
t0 = time.time()
print(f"[{args.quark} chunk {i}/{len(chunks)}] y = {chunks[i][0]:.4f} .. {chunks[i][-1]:.4f} ({len(chunks[i])} rows)", flush=True)
gen.run_command([str(ws.executable("fonllgridlha"))], d,
                gen.grid_input(prefix, pdf_id, var.mass, ffact=var.muF, fren=var.muR, y_values=chunks[i]),
                ws.logs / f"fonllgrid_{pdf_key}_{args.quark}_{var.tag}_n{args.nchunks}_chunk{i:02d}.log", ws.env())
if not chunk_complete(i):
    raise SystemExit(f"chunk {i} incomplete after fonllgridlha: {chunk_file(i)}")
print(f"[{args.quark} chunk {i}] done in {time.time() - t0:.0f} s", flush=True)
