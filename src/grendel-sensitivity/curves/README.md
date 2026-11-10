# Comparison curves and figures

Machine-readable reference curves (existing exclusions and the projections
of other experiments) for the HNL (BC6–8), dark-scalar (BC4 and BC5) and
fermiophilic-ALP (BC10) benchmarks, the renderer that draws GRENDEL's
sensitivity over them, and the digitisation tools that produced the curves.

```
data/hnl      HNLimits regions (transformed and raw), the single-Majorana BBN
              bound (cosmology/), projections, the BC7 compilation digitised
              two ways, per-experiment original sources, the HNLimits
              metadata workbook
data/bc4      PBC 2025 existing exclusion; SHiP (ECN3), CODEX-b projections
data/bc5      CODEX-b, MATHUSLA40, SHiP projections of the PBC 2025
              BC5 figure and ANUBIS of arXiv:2510.26932 (BR(h→SS) = 0.01; the PBC
              ANUBIS curve kept for reference); the exclusion is data/bc4's plus
              the ATLAS/CMS h→SS bands of arXiv:2509.02564 (existing/)
data/bc10     GKOZ LHCb recast, past beam dumps, NA62, FIPs 2022 and FIPs@LHCb 2025 bounds;
              SHiP, DarkQuest, CODEX-b, FASER2 projections
data/grendel_v1  the published GRENDEL tables of arXiv:2609.00152v1, a results
              tree (read-only); grendel_v2 those of release v2; the tables
              behind the figures of arXiv:2609.00152v2: grendel_v2_1_hybrid
              (BC4, BC5), grendel_v2_refined (HNL), grendel_v2_2 (BC10)
sources/      where the source figures come from (bundles not tracked)
plot/         the renderer
tools/        the extraction pipelines behind data/
```

## Rendering

```
pip install -e ".[curves]"      # matplotlib 3.10.9, HNLimits 1.2.0
python -m curves.plot paper --grendel-dir RESULTS --out figures
python -m curves.plot talk  --grendel-dir RESULTS --out figures [--with-envelope]
python -m curves.plot diagnostic --grendel-dir RESULTS --out figures
```

`RESULTS` (or `GRENDEL_RESULTS_DIR`) is the tree the scans write;
`curves/data/grendel_v1` holds the published tables in the same layout, so
`--grendel-dir curves/data/grendel_v1 --benchmark hnl bc4 bc10` renders the paper figures from a clean
checkout (`PROVENANCE.md` has their hashes; it has no `band/`):

```
RESULTS/hnl/sensitivity.csv    sensitivity_nsig10.csv    band/
RESULTS/bc4/sensitivity.csv    sensitivity_nsig10.csv    band/
RESULTS/bc5/sensitivity.csv    sensitivity_nsig10.csv
RESULTS/bc10/sensitivity.csv   sensitivity_nsig10.csv    band/
```

`sensitivity_nsig10.csv` is the `N_sig >= 10` island of the same scan
(`--thresholds 3 10`, split out by `python -m grendel.io.thresholds`); if
absent, the `N_sig >= 10` contour is omitted. `band/`
holds the variation-campaign outputs read by `--with-envelope` (BC4/BC10
single-source envelope) and by the HNL diagnostic figures.

`paper` writes the three manuscript figures at their final size with the
compact comparison set listed in `PROVENANCE.md`, plus the BC5 panel
(`bc5_grendel_paper`, not in the manuscript) in the same style; `talk`
uses the full diagnostic set (BC4/BC5/BC10; the HNL panels keep the paper
set) at projection size; `diagnostic` draws the HNL variation ribbons. The
HNL panels need the `HNLimits` package (its bundled data files are the
experimental bounds); BC4, BC5 and BC10 render without it.

## The curves

Every digitised `.dat` file carries its source, extraction method, date and
axis conventions in its header; the raw HNLimits copies carry their source,
mass unit and region type.
Three routes were used:

- **vector** (`tools/vector_*`): the contour is read from the figure's PDF
  path operators, calibrated from the axis ticks or frame. All BC4, BC5 and
  BC10 curves except the ALPINIST NA62 band, the BC7 compilation,
  CODEX-b's HNL panels and the BBN bound. (The arXiv source bundles carry
  the figures as vector PDFs and no data files, so this is the published
  curve itself, not a raster trace.)
- **raster** (`tools/raster_bc7`): an independent colour-mask digitisation
  of the BC7 compilation, kept as a cross-check on the vector route.
- **HNLimits** (`tools/hnlimits_pull`, `tools/hnlimits_process`): the
  per-experiment bounds and projections of mhostert/Heavy-Neutrino-Limits.
  `hnlimits_process` applies the workbook's unit, confidence-level and
  Dirac/Majorana rules and writes the regions the panels draw (the
  cosmology row replaced by the BBN bound);
  `hnlimits_process --check` regenerates them and compares bytes.

Policy: dedicated LLP experiments only, no ATLAS/CMS/LHCb projection
curves; existing exclusions are drawn as published (HNL bounds after
HNLimits' CL and Majorana normalisation). Conventions, the
per-panel layer choices and the handful of deliberate departures from
upstream data (single-Majorana BBN bound, FASER2 frame closure, CODEX-b
Majorana conversion, BC10 narrow-island guard, tick-geometry calibration)
are documented in `PROVENANCE.md`.

## Citing the sources

- PBC benchmark definitions: arXiv:1901.09966 (Beacham et al.)
- PBC BC5/BC7 figures: arXiv:2505.00947 (Alemany Fernández et al.)
- PBC ECN3 report, BC4 SHiP sensitivity: arXiv:2310.17726 (Ahdida et al.)
- HNLimits: github.com/mhostert/Heavy-Neutrino-Limits
- BBN bound on a single Majorana HNL: Chen and Zhang arXiv:2410.07343 (v1 drew
  Sabti et al. arXiv:2006.07387, which assumes two degenerate Majorana HNLs)
- SHiP HNL arXiv:1811.00930; FASER2 arXiv:1811.12522; CODEX-b arXiv:1911.00481;
  SHiP/DarkQuest ALP arXiv:2501.04525; GKOZ arXiv:2310.03524; FIPs 2022 arXiv:2305.01715
- ALPINIST (NA62 band): see `sources/alpinist/LICENSE` (BSD-3)
- per-experiment original sources: `data/hnl/original_sources/*/SOURCE.md`

The digitised files are derived from publicly available figures and data;
cite the originating papers when using them.
