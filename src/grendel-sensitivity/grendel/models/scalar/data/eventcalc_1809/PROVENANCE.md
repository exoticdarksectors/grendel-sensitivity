# EventCalc dark-scalar decay tables (Winkler arXiv:1809.01876 model), 0.01–63 GeV

`ctau.json.gz` is a list of `[m_S_GeV, ctau_m]` at `sin^2 theta = 1`;
`branching_ratios.json.gz` a list of `[channel, [pdg ids], [[m_S_GeV, BR], ...]]`
for 16 channels. They are the `scalar-1809` decay model of exHad
(V. Kryshtal, M. Ovchynnikov, arXiv:2609.16104; `data/scalar/eventcalc_{ctau,branching_ratios}_1809.01876.json`
at commit `7bb2308d1a60090603cbf07298d15b039f750727`; unchanged at `e41b07b`, the revision
pinned in `third_party/PINS.json`), which exHad in turn took from EventCalc / SensCalc
(M. Ovchynnikov et al., arXiv:2305.13383): the widths of arXiv:1809.01876
with the corrections and channels of arXiv:1904.10447. BSD-3-Clause; the
EventCalc notice is `LICENSE-EventCalc`.

Use in GRENDEL: BC5 needs the lifetime and the decay-mode composition up to
`m_h/2`, beyond the 7.5 GeV end of the digitised Winkler curves that fix the
published BC4 numbers (`../winkler_widths.csv`). `grendel.models.scalar.model`
keeps the digitisation below `WINKLER_M_MAX = 7.5 GeV` (BC4 unchanged) and
switches to these tables above it: `Gamma_tot = hbar c / ctau`, per-channel
widths `Gamma_tot * BR`, with `Jets-bb -> bb`, `Jets-cc -> cc`, `Jets-GG -> gg`,
`Jets-ss -> ss`, `tauPtauM -> tautau`, `muPmuM -> mumu`, `ePeM -> ee`. At the
seam (7.4 GeV) the two models differ by 16 % in `ctau` and by ~0.07 in
BR(tau tau); no smoothing is applied.

Copied 2026-09-25 from the exHad checkout in `grendel-mono/exhad/`.
