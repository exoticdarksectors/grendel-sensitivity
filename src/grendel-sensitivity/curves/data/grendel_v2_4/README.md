# GRENDEL v2.4 tables: BC10 with the corrected ALP widths

BC10 on the UP542 tunnel (geometry 2026-10-16.1) with the ALP-fermion widths of exHad e41b07b, whose
renormalization-group running from m_t to the GeV scale has the sign of the analytic NLL solution of arXiv:2012.12272
(the numerical Eq. (79) printed there has the opposite sign; the earlier exHad tables followed it). BC4, BC5 and the
HNL are unchanged: use `../grendel_v2_3` (exHad e41b07b changes only its ALP-fermion data). Do not edit these tables.

    python -m curves.plot paper --grendel-dir curves/data/grendel_v2_4 --benchmark bc10

| File | SHA-256 |
|---|---|
| `bc10/sensitivity.csv` | `5c08032e918fe1e0dcc7996c53374f16920cf061d6271295cb3c51d61dfc586d` |
| `bc10/sensitivity_nsig10.csv` | `db38e32f70338d6e61f0b7ffb124824be6743ee31c38c2cc8e4a791885db9c5e` |

## How they were made

* Templates: `grendel.models.alp.templates_exhad --n-templates 20000 --seed 1234 --workers 6 --sub-request 250
  --max-retries 10` with exHad e41b07b and Pythia 8.317 on macOS (exHad's 95 tests, including its stored event
  checksums, pass there), 151 masses (the eta and eta' windows are skipped), 11 sub-request retries, no failure.
  The mu mu convention factor is unchanged (16.518398).
* Four-vectors: those of `../grendel_v2_2/bc10` (GKOZ b -> s a coupling), plus 3.42-3.58 GeV in 20 MeV steps (no 3.50)
  with the same recipe (`grendel.models.alp.production --n-pool 1200000 --high-pt-tilt-scale 5 --seed 42`); controls at
  3.4 and 3.5 GeV reproduce the stored yields to 0.07 %.
* Scan: one HTCondor job per mass (153), one worker, the v2.3 code c621279, thresholds 3 and 10.

## Checks

* Against `../grendel_v2_3/bc10`: unchanged below 1.2 GeV (edges within 0.2 % at 0-0.6 GeV, median 1.002 at
  0.6-1.2 GeV); at 1.2-2 GeV lower edge median 1.023 (0.90-1.05), upper 1.27 (1.21-1.43); at 2-3 GeV lower 1.046,
  upper 1.21; at 3-3.6 GeV lower 0.98 (0.89-1.02), upper 1.32. Sensitive masses 117 (105), reach 3.5 GeV (3.34), the
  deepest point stays at 1.24 GeV (1/f = 1.036e-8 GeV^-1, 1.034e-8 before). The lifetime grows by x1.25-1.33 above
  1 GeV; peak N follows it exactly.
* Lifetime only (the v2.3 templates with the e41b07b lifetime, 18 masses): within 0.6 % of these rows in both edges up
  to 3.4 GeV (1.5-2 % at the marginal 3.5 GeV point), so the final-state change matters little.
* Control: the v2.3 templates rescanned with the same code reproduce the v2.3 rows bit for bit.
