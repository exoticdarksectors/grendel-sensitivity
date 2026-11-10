# GRENDEL HNL tables refined near the island tips

The published HNL tables (`../grendel_v1/hnl`, `../grendel_v2/hnl`) step by
0.2 GeV between 3.0 and 3.6 GeV and carry ~5 % Monte Carlo noise per mass on the
lower edge, which made the contours zigzag near the tips once the grid was
refined. These tables

- add 3.05, 3.10, 3.15, 3.25, 3.30, 3.35, 3.45, 3.50 and 3.55 GeV, and
- recompute every mass of 2.6-3.7 GeV (22 masses, Ue, Umu, Utau) with 4x
  production statistics;

all other rows are the published lines, verbatim. `grendel_v2_refined` (exHad
decays) draws the paper's HNL figure:
`python -m curves.plot paper --benchmark hnl --grendel-dir curves/data/grendel_v2_refined`.

| File | SHA-256 |
|---|---|
| `grendel_v1_refined/hnl/sensitivity.csv` | `e328b7dfe5bc6647c0dfe89d0684506c989969f520c4b14f0d08551328f66d46` |
| `grendel_v1_refined/hnl/sensitivity_nsig10.csv` | `7ade3638db3fb7985b04154774a2fef9a5f046d829ad7cf917f2163d35fb2b5a` |
| `grendel_v2_refined/hnl/sensitivity.csv` | `e52ac64b2c99485aa85c8bb146daa19085a67d6a66ea1eaf1bbaf321dfc374bb` |
| `grendel_v2_refined/hnl/sensitivity_nsig10.csv` | `c6ef30f88e4b7a5980f2b8b8f24edbd4df4025bc2faa188efc038c3d97500823` |

## The 2.6-3.7 GeV window (2026-10-06/07, `recipe/histat/`)

Production (this repository at d717b24, the `hnl` conda env of
`grendel-mono/hnl/environment.yml` for MadGraph + LHAPDF):

```
python -m grendel.models.hnl.production.run_all --flavor Ue Umu Utau --masses $MW \
    --channels bottom bc baryon wz --no-prompt-tau --skip-combine \
    --n-pool 400000 --wz-nevents 400000 --wz-nb-core 4 --workers 6 --seed 4242
```

then the closed channels (D mesons, kaons, induced tau) written empty, as in the
published production, and `combine_for_point`. v2: the vectors as produced
(SIGMA_BC_PB = 1.64e6 pb), exHad templates, scanned at d717b24. v1: B_c rows at
the v1 rate (x 0.9e6/1.64e6), FairShip templates, scanned with the v1.0 code
(b6670a8). Templates: the published sets where they have the mass, else those of
the refinement below. `recipe/histat/merge.py` replaces the window rows.

Why 4x: at 3.2-3.45 GeV the lower edge moves by 0.9 % rms with another scan seed,
1.4 % with another template set, but by 4.8 % rms (implied for the full sample)
between the two halves of the event sample: the noise is the production sample's
after reconstruction. With 4x statistics the lower-edge residuals around a smooth
trend fall from 6-12 % to 4-6 %, and what remains follows the Lambda_b and
B -> D l N thresholds (3.33, 3.41 GeV). The window rows moved by 7-8 % rms against
the 1x tables, with no systematic shift (mean +1 %).

## The added masses (2026-10-06, `recipe/`)

Production as above at 1x statistics (seed 42, `--n-pool 100000`, 100k W/Z
events). Templates for the new masses: FairShip
(`grendel.models.hnl.templates.fairship --n-templates 20000 --seed 1234`) in a
conda env with ROOT 6.38.04 and Pythia 8.312 (Homebrew's ROOT with Pythia 8.318
segfaults in TPythia8), and exHad (`grendel.models.hnl.templates.exhad
--n-templates 20000 --seed 1234`, exHad 7bb2308, Pythia 8.317). Checks:

- templates at 3.2 GeV regenerated with the same commands: FairShip byte-identical
  to the published v1 set for all three flavours; exHad identical in every array
  the scan reads;
- 3.2 and 3.4 GeV re-scanned from the published inputs: v2 reproduces the
  published rows exactly; v1 exactly except Ue 3.4 GeV (0.5 %), a row of the
  published v1 splice.

These masses are now inside the 4x window, so their 1x rows are superseded.
