# Historical ALPINIST/GKOZ inputs and production-RG provenance

The `digitized_width_<channel>.txt` files are the 2023 GKOZ decay-width
tables shipped by ALPINIST: Dalla Valle Garcia, Kahlhoefer, Ovchynnikov, and
Zaporozhchenko, arXiv:2310.03524, from
`widths/integrated_fermion_2310.03524/` under the BSD 3-Clause license.
`COMMIT_SHA.txt` pins the ALPINIST source revision fetched on 2026-07-08.

These digitized widths are retained for historical audit. They are **not** the
canonical decay input read by `model.py`. The current central model reads the
exact arXiv:2501.04525 SensCalc exports under `data/senscalc_2501/`; the named
2023 structural comparison reads the separately exported exact tables,
branching functions, and matrix elements under `data/senscalc_2310/`.

The historical table columns are `m_a [GeV]` and
`Gamma / (1/f)^2 [GeV^3]` in the BNT convention
`g_aff = c_f m_f/f`, with `1/f_BNT = 2/f_GKOZ`. Checks made when these files
were first adopted included `BR(a -> mu mu) ~ 9%` at 1 GeV and the charm-onset
structure near 2.58 GeV. Those checks describe the retired 2023 digitized
baseline; they are not validation targets for the 2025 central branching
fractions.

The directory remains relevant to production provenance. The effective
flavour-changing coefficient

```text
g_bs = CBS_EFF * (1/f),  CBS_EFF = 3.518383e-4
```

is evaluated for universal `c_f = 1` at `Lambda_UV = 1 TeV` with
`tools/compute_cbs_alpinist.py`, a BSD-attributed port of ALPINIST's
above-electroweak RG implementation of arXiv:2310.03524. There is no current
`CLL_RG` model constant; lepton and hadron decay information comes from the
pinned SensCalc exports described above.
