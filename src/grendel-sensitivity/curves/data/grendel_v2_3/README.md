# GRENDEL v2.3 tables: the grendel-geometry tunnel

The four benchmark tables on the UP542 tunnel of grendel-geometry v1.0 (geometry version 2026-10-16.1;
`GEOMETRY_UP542.md`), computed at c621279 (release-v2.3: 461b7e5 + the MadGraph helper fixes + the geometry mirror).
They are the record of these runs: do not edit them.

    python -m curves.plot paper --grendel-dir curves/data/grendel_v2_3 --benchmark bc4 bc5 bc10
    python -m curves.plot paper --grendel-dir curves/data/grendel_v2_3 --benchmark hnl      # needs HNLimits==1.2.0

| File | SHA-256 |
|---|---|
| `bc4/sensitivity.csv` | `ee52335f22c99e81d144540cb2214c9e6ce29d5f95074cdb8fdc02d8d1d0f80d` |
| `bc4/sensitivity_nsig10.csv` | `6414c78998f6849062addd6c89b42d849af6a81e979ae4d72bfe8c2d3cccaaa8` |
| `bc5/sensitivity.csv` | `e69e0882ea128b85b896201bd604c6d4c6cca8fa4a325f6120ed06619e8c0b47` |
| `bc5/sensitivity_nsig10.csv` | `f5318893dd3a4146651d49c60a9c8be6c40a587c21acdafadc0036d32dea3937` |
| `bc5/run_metadata.json` | `e34aabf32d610265ac298fef161a09d46db16043c0e702828e3a5f1e74fde3e1` |
| `bc5/paper_numbers.txt` | `c7e4b89fb757569d92d7b488e4c19ef2e9617ff2fdca8ed162bb878a139d2b2e` |
| `hnl/sensitivity.csv` | `58ceddf840126f706586f089f7bc76d4898c281dc5f8853d24f66e7957644082` |
| `hnl/sensitivity_nsig10.csv` | `2f2dc81011c56d5b96a8edf2dfd6d70839ee640d07f0110d3520e37a96a13a9f` |
| `bc10/sensitivity.csv` | `bcd23140c709426064b9feba2c57b843916639816f0eb45278a861a073d3a8fc` |
| `bc10/sensitivity_nsig10.csv` | `8c6d43e642fddcd2d73c91ed8fad14e4871bcd3e2495612c68297295002b6b64` |

## How they were made

One HTCondor job per (benchmark, flavour, mass), one worker each, with the scan commands of the earlier tables; a
scan gives every mass its own seed, so a one-mass job writes the row a multi-worker scan writes. The four-vectors and
templates are those of the earlier tables (they do not depend on the tunnel):

* **BC4**: the 86 masses of `../grendel_v1/bc4`, the v2 four-vector pool and the hybrid templates of
  `../grendel_v2_1_hybrid` (Winkler lifetime and split, exHad final states), 100 decay samples, thresholds 3 and 10.
* **BC5**: the 112 masses of `../grendel_v2_1_hybrid/bc5`, `scan_bc5 --decay-samples 100`, the same hybrid templates;
  four-vectors with the Higgs (pT, y) histogram showered with Pythia 8.317 (461b7e5).
* **HNL**: production with HNLCalc b1d36c6 (FLAG 2024 decay constants): charm, kaon (Pythia 8.317 spectrum) and B_c
  regenerated (seed 42), the other channels of the v2 production; exHad templates of v2. 0.2-2.0 and 4.2-10 GeV at
  the v2 masses (104 per flavour); 2.0-4.2 GeV a fine grid of 137 masses per flavour (25 MeV steps, 5 MeV at
  3.6-3.9 GeV) at 10x the published statistics: B, B_c, Lambda_b pools of 1 000 000 (`run_all --channels bottom bc
  baryon`, seeds 7000/7100/7200 for Ue/Umu/Utau), 1 000 000 MadGraph W/Z events per point, exHad templates of
  20 000 decays (seed 1234). Every 2.0-4.2 GeV row of the coarse grid is replaced.
* **BC10**: the v2.2 inputs of `../grendel_v2_2/bc10` (GKOZ b -> s a coupling, exHad templates), the default decay
  samples. Interim: the ALP widths are being revised upstream (exHad and SensCalc), and BC10 follows them.

## Checks

* Against the same scans on geometry 2026-10-14.2 (identical except sub-mm in the centreline and a few mm in the
  section): medians per mass range within 0.4 % for every benchmark; single masses BC4 0.996-1.005, BC5 lower edge
  0.994-1.010 and upper 0.986-1.027, HNL 0.986-1.034 at 0.5-3.4 GeV and 0.95-1.03 below 0.5 GeV and at the island
  tips (statistics); reach and sensitive masses unchanged (`comparison/summary.md` of the campaign).
* Linux (lxplus) against macOS for the same code and inputs (36 BC4 masses, an HNL point): every edge column bit for
  bit; peak N and the effective-sample-size diagnostics agree to 1e-14 (summation order).
* The scan-only pipeline with the code of 461b7e5 on the stored four-vectors reproduces the 14 Oct row of Umu 3 GeV
  byte for byte.
* Against the v2.2 tables (the arXiv:2609.00152 tunnel, 1.5-3.1 m too low and mirrored): hits x0.89-0.91; lower edge
  x1.05 (BC4), x1.04-1.11 (HNL below 3.4 GeV), x1.07-1.11 (BC5), x1.02-1.06 (BC10); upper edge x0.83-0.94; BC4 reach
  3.8 -> 3.7 GeV; HNL closure Ue 3.66, Umu 3.665, Utau 3.855 GeV (v2.2: 3.7, 3.7, 4.0 on its coarser grid).
* `bc5/paper_numbers.txt`: the quantities Sec. IV.C quotes.
