# GRENDEL v1 contour tables

The six frozen GRENDEL tables drawn in Figs. 4–6 of arXiv:2609.00152v1,
copied unchanged (read-only) from the analysis' publication freeze. They
are a results tree in the layout the renderer reads, so the figures render
from a clean checkout:

```
python -m curves.plot paper --grendel-dir curves/data/grendel_v1 --benchmark hnl bc4 bc10 --out figures
```

| File | Figure | Published as | SHA-256 |
|---|---|---|---|
| `bc4/sensitivity.csv` | Fig. 4, `N_sig >= 3` | `bc4_island.csv` | `45aec445068cfe77ec2d75bf20e251255559545561cde3f77f1c2ba74c35e84c` |
| `bc4/sensitivity_nsig10.csv` | Fig. 4, `N_sig >= 10` | `bc4_island_nsig10.csv` | `515c04302f26a2241123a400e9b88566e32dd2f5f587a11fa7f52a6c0495c617` |
| `hnl/sensitivity.csv` | Fig. 5, `N_sig >= 3` | `grendel_hnl_sensitivity.csv` | `b051ac77b1c2dd677c809e6acd6f3a34c8db374fb0d474a775d1678777a9b407` |
| `hnl/sensitivity_nsig10.csv` | Fig. 5, `N_sig >= 10` | `grendel_hnl_sensitivity_nsig10.csv` | `ca50461b92b4ecc8586de2faacac3abf6248e9ba3647556e32e31249b22f2b93` |
| `bc10/sensitivity.csv` | Fig. 6, `N_sig >= 3` | `bc10_sensitivity.csv` | `d4ea435258a4cdecef47b5764f36bff452d2f08f6e3190172a3f1e97f7ec308c` |
| `bc10/sensitivity_nsig10.csv` | Fig. 6, `N_sig >= 10` | `bc10_sensitivity_nsig10.csv` | `789aecb1bedc611db2e5ac16a51813c5e18fb633088200b3bca916dec02661c7` |

They are the published record, not outputs this repository regenerates:
do not edit them, and point new scans at their own results tree. The HNL
`N_sig >= 3` table is a splice of two scan campaigns (the `K -> l N`
window, `Ue m_N <= 0.485` and `Umu m_N <= 0.38` GeV, from the later one);
re-solved in one campaign it differs on 40 of 369 rows, by at most 2.7 %.
Its `N_sig >= 10` partner is that self-consistent re-run, nested inside it.

With the v1 comparison layers these tables reproduce the published PNGs
byte for byte; `curves/PROVENANCE.md` lists both sets of hashes and which
layers have changed since.
