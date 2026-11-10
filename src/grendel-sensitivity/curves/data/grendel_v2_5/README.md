# GRENDEL v2.5 tables: BC4 and BC5 with the EventCalc scalar model

BC4 and BC5 on the UP542 tunnel (geometry 2026-10-16.1) with one scalar decay model at every mass: the lifetime,
branching fractions and final states of exHad's `scalar-1809` model, i.e. the widths of arXiv:1809.01876 with the
corrections and channels of arXiv:1904.10447 as tabulated in EventCalc/SensCalc (arXiv:2305.13383; the tables of
`grendel/models/scalar/data/eventcalc_1809`). It is the model of the SHiP BC4 and BC5 projections. `../grendel_v2_3`
used the Winkler digitisation below 7.5 GeV (this table above it). BC10 is `../grendel_v2_4`, the HNL
`../grendel_v2_3`. Do not edit these tables.

    python -m curves.plot paper --grendel-dir curves/data/grendel_v2_5 --benchmark bc4 bc5

| File | SHA-256 |
|---|---|
| `bc4/sensitivity.csv` | `d0dcf097bbfab2a2291d0a081959e69a49143c3e5807d9b2ab46a41cd62f77d4` |
| `bc4/sensitivity_nsig10.csv` | `08aa20d673965aa2c6b57d3a946bb6b349e466b7d7c3270d6e98b4d9256f9197` |
| `bc5/sensitivity.csv` | `243a529a117953249b694e283dd2eac3161661640137a27b8f73991f2c62a2df` |
| `bc5/sensitivity_nsig10.csv` | `12ce40782150df701f68936acaafd794bc9b4c22f21444c635c53503a9125123` |
| `bc5/run_metadata.json` | `984df8f9c8047e7932a8451b57dc462166959df72a436b1307f47cae5dfbea77` |
| `bc5/paper_numbers.txt` | `72b97db25d894230da4f729e85f58096520f471aa169ece8c7c7006169b9ad28` |

## How they were made

Exactly the runs of `../grendel_v2_3` (code c621279, the same four-vectors, masses, seeds, 100 decay samples for
BC5, thresholds 3 and 10, one HTCondor job per mass), with only the decay templates changed: the exHad bundles from
which `../grendel_v2_1_hybrid/make_hybrid_templates.py` built the hybrid ones (exHad 7bb2308, whose scalar data are
unchanged at the pinned e41b07b; Pythia 8.317; 20 000 per mass, seed 1234), used as they are, so `ctau_m_u2eq1` and
the channel mix are exHad's. `bc5/paper_numbers.txt` is
`grendel/campaigns/bc5/bc5_v2_1_20261003/results/bc5_hybrid/paper_numbers.py` on this table.

## Checks

* Lifetime, EventCalc / Winkler: 1.03 below 1 GeV, except near the f0(980), where the lifetime minimum lies at
  0.975 GeV in Winkler's model and at 1.0 GeV in EventCalc's (ratio 1.39 at 0.975 GeV, 0.68 at 1.0 GeV); 0.83 at
  1-1.5 GeV, 0.59 at 1.5-2 GeV, 0.48 at 2-3.5 GeV, 0.44 above 3.5 GeV in BC4's range, 1.07 at 5-10 GeV and 1.00
  above 10 GeV.
* BC5 against `../grendel_v2_3/bc5`: both edges follow the lifetime ratio (Higgs and quartic B production do not
  depend on the mixing), medians 1.02, 0.84, 0.60, 0.50, 0.88, 1.07 and 1.00 in those ranges; 112/112 masses
  sensitive at N >= 3 and N >= 10; deepest point 2.05e-17 at 50 GeV (2.06e-17).
* BC4 against `../grendel_v2_3/bc4`: lower edge 1.01 below 1 GeV, 0.92 at 1-1.5 GeV, 0.79 at 1.5-3.5 GeV; upper
  edge 0.78, 0.49, 0.36 and 0.24 in the ranges above 1 GeV. 75 masses sensitive (76), reach 3.6 GeV (3.7), N >= 10
  to 3.0 GeV (3.6); deepest point 7.2e-12 at 1.0 GeV (7.0e-12 at 0.975 GeV), where the EventCalc lifetime has its
  f0(980) minimum.
* Which width each comparison curve uses (`grendel/campaigns/bc10/bc10_v24_20261019/codexb_check`): the SHiP BC4
  (Briefing Book) and BC5 (PBC 2025) upper edges are smooth in c tau only with the EventCalc lifetime
  (`upper_edge_lifetime_test.*`); the PBC 2025 BC5 curves of CODEX-b, MATHUSLA40 and ANUBIS imply lifetimes 3.0-3.3,
  1.6-1.7 and 1.1-1.2 times shorter than this one at 15 GeV (`implied_lifetime_bc5.*`).
