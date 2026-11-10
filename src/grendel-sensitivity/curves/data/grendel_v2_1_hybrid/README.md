# GRENDEL scalar tables: exHad final states with the Winkler lifetime

The BC4 and BC5 tables behind Figs. 4 and 5 of the v2 paper (2026-10-07). The
decays keep the Winkler lifetime and the split between e e, mu mu, tau tau and
hadronic decays of v1 (the PBC convention, and that of every comparison curve)
and take the final states from exHad. They are the record of these runs: do not
edit them.

    python -m curves.plot paper --grendel-dir curves/data/grendel_v2_1_hybrid --benchmark bc4 bc5

| File | SHA-256 |
|---|---|
| `bc4/sensitivity.csv` | `c6152e485b2fe71c6a1727c96c9d6867a808497ab5b5ae8dd858fa6f0678f05e` |
| `bc4/sensitivity_nsig10.csv` | `61794491dd1ba2fff1ead3e263718e2063ea4585d92023e887a5865e0ad36bc3` |
| `bc5/sensitivity.csv` | `4a06fb915d5cd3c665bbc028afee161afc6e354f105d4b98579b900bca71daf8` |
| `bc5/sensitivity_nsig10.csv` | `61a4478199652b3f54e9620a99bfc9cdf759bbe0766fd7f309708af51dc569f3` |
| `bc5/run_metadata.json` | `34f74d8e161b88efe8f3c1a4e9a642da538a6d4e33987e6706602e5a56fe056e` |
| `bc5/paper_numbers.txt` | `2e535a99f1382b57cc883773a21b4b601efdbff38b9ff7db741cf0540c023c2f` |

## Decay templates

`make_hybrid_templates.py` (SHA-256 `64c65460ec82eec3e1884ef7cd4c0e246358974cf9d12a9c48779c8aa2ca4c1a`) builds them from the exHad
`scalar-1809` bundles of the v2 scalar templates
(`grendel/campaigns/v2_baseline_20260929/templates/scalar`, exHad 7bb2308,
Pythia 8.317, 20 000 per mass): `ctau_m_u2eq1` is set to the Winkler value the
bundle stores as `ctau_m_u2eq1_reference` (equal to
`model.ctau_sin2theta1(m, 'winkler')`, which above 7.5 GeV continues with the
EventCalc tables of arXiv:1809.01876 with the arXiv:1904.10447 corrections), and
the events are subsampled without replacement (seed 20261007) so that the e e,
mu mu, tau tau and hadronic fractions are the Winkler branching ratios; inside
the hadronic group the composition and every final state are exHad's. 8 850 to
20 000 templates per mass (median 18 064).

    python make_hybrid_templates.py <v2 scalar templates> grendel_work/bc5/templates_hybrid

## BC5 (grendel_light clone at 11fd035 = v2.1)

    GRENDEL_WORK_DIR=grendel_work python -m grendel.models.scalar.scan_bc5 --workers 5 \
      --decay-samples 100 --thresholds 3 10 --vectors-dir grendel_work/bc5/llp_4vectors \
      --templates-dir grendel_work/bc5/templates_hybrid --geometry-dir results/bc5/geometry_cache \
      --out results/bc5_hybrid --mass <the interleaved grid of ../grendel_v2_1/README.md>
    python -m grendel.io.thresholds results/bc5_hybrid/sensitivity.csv --threshold 10

The v2.1 run of `../grendel_v2_1/bc5` (same code, four-vectors, geometry cache, per-mass seeds,
100 decay samples, thresholds 3 and 10) with only the decay stage changed. Run in two parts
(35 + 77 masses, `--resume`), 112/112 sensitive at N >= 3 and N >= 10.

### Checks

* Hit counts and lifetimes equal those of `../grendel_v2_1/bc5` at every mass: only the decay stage differs.
* Above 8 GeV, where the Winkler and exHad lifetimes coincide, the edges equal the full-exHad table
  `../grendel_v2/bc5` to a median 0.006 dex (lower) and 0.002 dex (upper).
* Against `../grendel_v2_1/bc5` (two-body proxy): identical below 0.6 GeV; lower edge -0.13 dex at
  1.3-2 GeV (four-pion and fragmentation final states are more often visible than the proxy's),
  +0.05 to +0.08 dex at 2-30 GeV (the proxy counted every continuum decay as two charged tracks),
  -0.09 dex at 30-63 GeV; upper edge within 0.09 dex everywhere.
* `bc5/paper_numbers.txt`: the quantities Sec. IV.C quotes (`results/bc5_hybrid/paper_numbers.py`).

## BC4 (this repository at b7aad98)

The v2 BC4 inputs (`../grendel_v2/README.md`: the v1 four-vector pool, `--n-pool 400000
--high-pt-tilt-scale 5`, seed 42; the v2 geometry cache) on the 86 masses of the v1
table, with the BC5 template set above (it covers every BC4 mass):

    python -m grendel.models.scalar.scan --workers 4 --resume --thresholds 3 10 \
      --vectors-dir <BC4 pool> --templates-dir grendel_work/bc5/templates_hybrid \
      --geometry-dir <v2 BC4 geometry cache> --out RESULTS/bc4 \
      --mass <the masses of ../grendel_v1/bc4/sensitivity.csv>
    python -m grendel.io.thresholds RESULTS/bc4/sensitivity.csv --threshold 10

100 decay samples, p > 100 MeV, the per-mass seeds (the scan defaults); 77 of 86 masses
sensitive. The appendix cutflow column (0.975 GeV, the deepest reach):

    python -m grendel.cutflow --model bc4 --mass 0.975 --from-results curves/data/grendel_v2_1_hybrid/bc4/sensitivity.csv \
      --vectors-dir <BC4 pool> --templates-dir grendel_work/bc5/templates_hybrid --geometry-dir <geometry cache>

### Checks

* The same scan without `--templates-dir` (the analytic decays of v1), run at 0.975 and 2 GeV, reproduces
  `../grendel_v1/bc4` to 0.5 %, and its cutflow at the v1 deepest reach the v1 appendix
  column to 0.001.
* Against `../grendel_v1/bc4`: identical below 1 GeV (deepest 6.50e-12 at 0.975 GeV); lower
  edge x0.82-0.86 at 1.1-1.7 GeV, where four-pion final states are more often visible than
  v1's charged-track models (acceptance x1.45, as BC5's -0.13 dex there); both edges 1-16 %
  lower at 2-3.8 GeV; reach 3.8 GeV (N >= 3) and 3.6 GeV (N >= 10), unchanged.
