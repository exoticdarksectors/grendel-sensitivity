# Comparison curves and figures

Machine-readable reference curves (existing exclusions and the projections
of other experiments) for the HNL (BC6–8), dark-scalar (BC4) and
fermiophilic-ALP (BC10) benchmarks, the renderer that draws GRENDEL's
sensitivity over them, and the digitisation tools that produced the curves.

```
data/hnl      HNLimits regions (transformed and raw), projections, the BC7
              compilation digitised two ways, per-experiment original sources,
              the HNLimits metadata workbook
data/bc4      PBC 2025 existing exclusion; SHiP, CODEX-b projections
data/bc10     GKOZ LHCb recast, past beam dumps, NA62, FIPs 2022 bounds;
              SHiP, DarkQuest projections
sources/      where the source figures come from (not tracked)
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

`RESULTS` (or `GRENDEL_RESULTS_DIR`) is the tree the scans write:

```
RESULTS/hnl/sensitivity.csv    sensitivity_nsig10.csv    band/
RESULTS/bc4/sensitivity.csv    sensitivity_nsig10.csv    band/
RESULTS/bc10/sensitivity.csv   sensitivity_nsig10.csv    band/
```

`sensitivity_nsig10.csv` is the same scan solved at `N_sig >= 10`
(`--thresholds 3 10`); if absent, the dashed contour is omitted. `band/`
holds the variation-campaign outputs read by `--with-envelope` (BC4/BC10
single-source envelope) and by the HNL diagnostic figures.

`paper` writes the three manuscript figures at their final size with the
compact comparison set listed in `PROVENANCE.md`; `talk` uses the full
diagnostic set at projection size; `diagnostic` draws the HNL variation
ribbons. The HNL panels need the `HNLimits` package (its bundled data
files are the experimental bounds); BC4 and BC10 render without it.

## The curves

Every `.dat` file carries its source, extraction method, date and axis
conventions in its header. Three routes were used:

- **vector** (`tools/vector_*`): the contour is read from the figure's PDF
  path operators, calibrated from the axis ticks. All BC4 and BC10 curves,
  the BC7 compilation, and CODEX-b's HNL panels.
- **raster** (`tools/raster_bc7`): an independent colour-mask digitisation
  of the BC7 compilation, kept as a cross-check on the vector route.
- **HNLimits** (`tools/hnlimits_pull`, `tools/hnlimits_process`): the
  per-experiment bounds and projections of mhostert/Heavy-Neutrino-Limits.
  `hnlimits_process` applies the workbook's unit, confidence-level and
  Dirac/Majorana rules and writes the regions the panels draw;
  `hnlimits_process --check` regenerates them and compares bytes.

Policy: dedicated LLP experiments only, no ATLAS/CMS/LHCb projection
curves; existing exclusions are drawn as published. Conventions, the
per-panel layer choices and the handful of deliberate departures from
upstream data (per-flavour cosmology, FASER2 frame closure, BC10
narrow-island guard, tick-geometry calibration) are documented in
`PROVENANCE.md`.

## Citing the sources

- PBC benchmark definitions: arXiv:1901.09966 (Beacham et al.)
- PBC BC5/BC7 figures: arXiv:2505.00947 (Antel et al.)
- PBC ECN3 report, BC4/BC10 sensitivities: arXiv:2310.17726 (Ahdida et al.)
- HNLimits: github.com/mhostert/Heavy-Neutrino-Limits; cosmology Sabti et al. arXiv:2006.07387
- SHiP HNL arXiv:1811.00930; FASER2 arXiv:1811.12522; CODEX-b arXiv:1911.00481;
  SHiP/DarkQuest ALP arXiv:2501.04525; GKOZ arXiv:2310.03524; FIPs 2022 arXiv:2305.01715
- ALPINIST (NA62 band): see `sources/alpinist/LICENSE` (BSD-3)
- per-experiment original sources: `data/hnl/original_sources/*/SOURCE.md`

The digitised files are derived from publicly available figures and data;
cite the originating papers when using them.
