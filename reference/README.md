# Reference material

The validation steps compare what the pipelines rebuild with what the original runs produced. No production step
reads any of it: without it the tables and figures come out the same, and fewer comparisons are made.

| What | Where | Used for |
|---|---|---|
| The published tables | `src/grendel-sensitivity/curves/data/` (`$GRENDEL_PUBLISHED`) | every rebuilt table is compared row by row |
| The paper's figures | `reference/paper_figures/` (`$GRENDEL_PAPER_FIGURES`): the PNGs of Figs. 4-7 and the PDFs of Figs. 4, 5 and 7 | pixel and byte comparisons of the figures |
| The outputs of the original runs | not in this repository; `$GRENDEL_REF` | everything else: four-vector pools, MadGraph LHE files, exHad template bundles, ray-casting caches, SensCalc runs, FONLL workspace, campaign records |

## The outputs of the original runs (`GRENDEL_REF`)

`files/<pipeline>.txt` lists every file of the original runs that a validation step reads (paths relative to the
root of the reference tree). The lists were recorded by running every validation step under `strace`.
`MANIFEST.tsv` gives the size and SHA-256 of each:

| Pipeline | Files | Size | Mostly |
|---|---|---|---|
| `scalar` | 241 | 16.2 GB | the stored BC5 four-vector pool (15 GB), the published exHad bundles, the 2026-09-25 MadGraph LHE files and their shower, the BC4 pool manifest |
| `hnl` | 1974 | 18.1 GB | the as-built coarse-grid four-vectors (1248 files), the fine-grid MadGraph records (412), the lxplus scan rows (312), the template manifests |
| `alp` | 467 | 1.05 GB | the exHad bundles of the v2.4 and v2.0 runs (294), the ray-casting caches (18), the run records |
| `senscalc` | 33 | 0.67 GB | the SensCalc eb2c456 tarball the campaign used (0.6 GB), the campaign's SensCalc outputs and logs, M. Ovchynnikov's SHiP file, the July 2026 export |
| `fonll` | 11 | 7 MB | the grids and calibration of the original FONLL workspace |
| all | 2726 | 36.0 GB | |

488 files (CSV files of the original runs, 487 of them HNL four-vectors) are empty in the reference tree the
reproduction was made with: they lost their content when copied from an HFS+-compressed macOS volume. The
validations compare those against the sizes and digests the original runs recorded, not against their content.

**Using it.** Put the tree anywhere and point `GRENDEL_REF` at its root before `source env/activate.sh`:

```
export GRENDEL_REF=/path/to/reference-tree
source env/activate.sh
python tools/reference_manifest.py --check      # every listed file present, with the size and SHA-256 above
```

Without `GRENDEL_REF` (or with `GRENDEL_REF=none`) every comparison that needs it prints a `SKIP` line, also
collected in `work/skipped_validations.txt`, and the pipeline carries on. On the machine the reproduction was made
on, the tree is `../extra_space`, the default when it exists.

**Getting it.** `tools/make_reference_bundle.sh DEST` packs the listed files into one archive per pipeline, which
unpack into a reference tree. Those archives are meant for a public data record (e.g. on Zenodo); until one
exists, ask the authors. `not_redistributable.txt` names the files the archives leave out because they were given
to the authors privately: M. Ovchynnikov's SensCalc result for SHiP (2026-10-16), from which the published SHiP curve
of Fig. 7 was converted. Without it, `pipelines/senscalc` uses the published curve and skips that comparison.
