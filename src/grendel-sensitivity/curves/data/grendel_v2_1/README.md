# GRENDEL v2.1 BC5 table

The BC5 table behind the BC5 figure of the v2 paper draft (the island split
by production mode and rescaled to BR(h -> SS) = 0.1 %). It is the record of
that run: do not edit it.

| File | SHA-256 |
|---|---|
| `bc5/sensitivity.csv` | `157c3c23256f16454ada28f4670d9949a40ce016cf07ef7f025c7ea777a1c1ea` |
| `bc5/sensitivity_nsig10.csv` | `398d46b32ad5236554d48309c1e5547d23ec056ab40dd98a7b5e070f3d4c515e` |
| `bc5/run_metadata.json` | `4636428587d025b166cf6a79da4d1c960c23a43e6d4ddc63f06471285e399055` |
| `bc5/paper_numbers.txt` | `7a0468e4661fb9e380b6c34dfe8e1f64296a6eb6d26a8f94a48b4cf2928a95b2` |

`run_metadata.json` carries the SHA-256 of every input four-vector file;
`paper_numbers.txt` the numbers the paper text quotes. The per-mode tables
are not stored: `python -m grendel.io.variants bc5/sensitivity.csv --variant
<hSS|mixing|BSS|brhss0.001|brhss0.001_hSS|brhss0.001_BSS>` and
`python -m grendel.io.thresholds bc5/sensitivity.csv --threshold 10`
regenerate them, and the N >= 10 table, byte for byte.

## How it was made (2026-10-03)

Code: v2.1 (this repository from 8cd8077; public tag `v2.1`, 11fd035). The
run used 5a00487 and its four-vectors c82e40e, which differ from 11fd035 only
in `scripts/plot_bc5_modes.py`. Two deviations from the v2 BC5 table
(`../grendel_v2/bc5`): the analytic decay model (`winkler` widths) instead
of the exHad templates, and the masses scanned in the interleaved order
below. 100 decay samples, thresholds 3 and 10, scan_bc5 production defaults
(seed 42), `GRENDEL_WORK_DIR=grendel_work`:

```
M="0.14 0.18 0.22 0.26 0.3 0.34 0.38 0.42 0.46 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.2 2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4 4.2 4.4 4.6 4.8 5.5 7 9 11 14 16 20 25 30 40 50 60 0.16 0.2 0.24 0.28 0.32 0.36 0.4 0.44 0.48 0.525 0.575 0.625 0.675 0.725 0.775 0.825 0.875 0.925 0.975 1.05 1.15 1.25 1.35 1.45 1.55 1.65 1.75 1.85 1.95 2.1 2.3 2.5 2.7 2.9 3.1 3.3 3.5 3.7 3.9 4.1 4.3 4.5 4.7 5 6 8 10 12 15 18 22.5 27.5 35 45 55 62"
python -m grendel.models.scalar.scan_bc5 --workers 12 --decay-samples 50  --thresholds 3 10 --out results/bc5 --mass $M
python -m grendel.models.scalar.scan_bc5 --workers 12 --decay-samples 100 --thresholds 3 10 --out results/bc5_ds100 \
    --geometry-dir results/bc5/geometry_cache --mass $M
```

The first run writes the four-vectors and the geometry cache the second
reuses; this table is the second.

## The figure

```
python scripts/plot_bc5_modes.py curves/data/grendel_v2_1/bc5/sensitivity.csv --out bc5_modes.png
```

gives the paper's `bc5_grendel_paper.png` byte for byte with matplotlib
3.11.2, numpy 2.5.3, pandas 3.0.6 (Python 3.13); other matplotlib versions
differ by a pixel row.
