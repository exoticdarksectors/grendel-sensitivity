# GRENDEL v2.2 BC10 table: the b -> s a coupling of arXiv:2310.03524

The BC10 table behind Fig. 7 of the v2 paper (2026-10-05). It is the v2 BC10
run (`../grendel_v2/bc10`: exHad decays, the same templates, geometry, mass grid
and scan settings) with the production coupling of arXiv:2310.03524 Table 1
(`model.CBS_EFF` = 4.405e-4 instead of ALPINIST's 3.518e-4, CHANGELOG 2.2).
It is the record of that run: do not edit it.

    python -m curves.plot paper --benchmark bc10 --grendel-dir curves/data/grendel_v2_2

| File | SHA-256 |
|---|---|
| `bc10/sensitivity.csv` | `20eeda6b9ada348dfb4b4ac328a2e80b0b3b92a075bcc5d6821328e8d4790abf` |
| `bc10/sensitivity_nsig10.csv` | `5bede2a4f2f0f827bbdc5995173d1690b33954cda702e30a5cde44eaec5717d3` |

## How it was made

Code: this repository at 73f1203 (the coupling commit).

    python -m grendel.models.alp.production --n-pool 1200000 --high-pt-tilt-scale 5 --out-dir VECTORS
    python -m grendel.models.alp.scan --workers 6 --thresholds 3 10 --vectors-dir VECTORS \
        --templates-dir <v2 exHad BC10 templates> --geometry-dir <v2 BC10 geometry cache> \
        --out RESULTS --mass <the 145 masses of ../grendel_v2/bc10/sensitivity.csv>
    python -m grendel.io.thresholds RESULTS/bc10/sensitivity.csv --threshold 10

Production seed 42 (the default); 60 decay samples, p > 100 MeV, the per-mass
seeds (the scan defaults). Templates: exHad 7bb2308 with Pythia 8.317, 20 000 per
mass, base seed 1234, as `../grendel_v2/README.md`. 2353 s on 6 workers.

The appendix cutflow column (1.49 GeV, 1/f = 9.76e-9 GeV^-1, the deepest reach):

    python -m grendel.cutflow --model bc10 --mass 1.49 --from-results curves/data/grendel_v2_2/bc10/sensitivity.csv \
        --vectors-dir VECTORS --templates-dir <templates> --geometry-dir <geometry cache>

## Checks

* The four-vectors are those of v2 (the v1 production) with every weight x1.567236
  ((4.405/3.518)^2) and identical kinematics.
* Against `../grendel_v2/bc10`: peak N x1.5672 at every mass; lower edge x0.87
  (median; 0.63-0.89), upper edge x1.11 (1.03-1.53); sensitive masses 96 -> 112
  (gains at 1.25, 1.38-1.42, 1.47 and 3.24-3.40 GeV).
* The paper figure leaves out the 1.374-1.423 GeV island (49 MeV wide, peak N 3.97):
  the renderer drops islands narrower than 50 MeV.
