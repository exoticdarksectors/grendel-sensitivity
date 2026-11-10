# HNL validation (work/hnl/validation)

## 1. Coarse-grid channels (step 05, 06)

Dmeson       identical files 208/312; rows text-identical 15399668/15399866; rows with kinematics differing 198 (max rel diff 3.7e-08); w_sum ratio 1.0000000-1.0000000
Bc           identical files 163/312; rows text-identical 26199771/26200000; rows with kinematics differing 229 (max rel diff 1.6e-08); w_sum ratio 1.0000000-1.0000000
Kmeson       identical files 279/312; rows text-identical 9801211/9900000; rows with kinematics differing 98788 (max rel diff 9.4e-05); w_sum ratio 1.0000000-1.0000000
Bmeson       identical files 133/312; rows text-identical 24899440/24899755; rows with kinematics differing 315 (max rel diff 2.1e-08); w_sum ratio 1.0000000-1.0000000
Bbaryon      identical files 156/312; rows text-identical 22299755/22300000; rows with kinematics differing 245 (max rel diff 2.6e-08); w_sum ratio 1.0000000-1.0000000
induced_tau  identical files 154/312; rows text-identical 21099491/21099789; rows with kinematics differing 298 (max rel diff 1.7e-08); w_sum ratio 1.0000000-1.0000000
tau          211 non-empty points: identical files 0, rows text-identical 22835556/34384349; w_sum ratio median 1.000000 (1.000000-1.000000); <E> ratio median 0.99994
WZ           312 non-empty points: identical files 0, rows text-identical 21645202/31200000; w_sum ratio median 1.000000 (1.000000-1.000000); <E> ratio median 1.00000

## 2. Fine-grid production (scan jobs, step 06)

Bmeson: rows identical at 411/411 points; weight sum identical to every printed digit at 76, max |new/published - 1| 4.4e-11
Bc: rows identical at 411/411 points; weight sum identical to every printed digit at 0, max |new/published - 1| 2.8e-11
Bbaryon: rows identical at 411/411 points; weight sum identical to every printed digit at 307, max |new/published - 1| 2.6e-11
W/Z weight sum (sigma x 1.3) new/published at 411 points: median 1.00001, range 0.99995-1.00006
MadGraph cross section at 411 points: ratio median 1.00000, pull mean +0.04 rms 0.62

## 3. Templates (step 07)

696 templates compared
ctau(U^2=1) identical at 696/696; max rel diff 0.0e+00
channel fractions: chi2/ndf summed 0/16098; p-values: 0 below 0.01 (expected 7.0)
mean multiplicity ratio 1.0000 (0.9978-1.0022); visible fraction diff -0.0000 (rms 0.0001, binomial sigma ~0.0040)

## 4. Rows (step 08, 09)

723 of 723 rows; sensitivity flag (N>=3) agrees at 721, (N>=10) at 723
| flavour | mass range | rows | u2_min new/pub median (min-max) | u2_max | peak N |
|---|---|---|---|---|---|
| Ue | 0-0.5 | 20 | 1.000 (0.992-1.000) | - | 1.000 (1.000-1.001) |
| Ue | 0.5-2.0 | 60 | 1.000 (0.998-1.003) | 1.000 (0.994-1.002) | 1.000 (0.996-1.003) |
| Ue | 2.0-3.4 | 56 | 1.000 (0.997-1.006) | 0.999 (0.986-1.010) | 1.000 (0.988-1.015) |
| Ue | 3.4-4.2 | 80 | 1.001 (0.881-1.206) | 0.996 (0.860-1.072) | 0.999 (0.973-1.024) |
| Ue | 4.2-11 | 25 | - | - | 1.000 (0.985-1.010) |
| Umu | 0-0.5 | 20 | 1.000 (1.000-1.000) | - | 1.000 (1.000-1.000) |
| Umu | 0.5-2.0 | 60 | 1.000 (0.997-1.001) | 1.000 (0.999-1.002) | 1.000 (1.000-1.007) |
| Umu | 2.0-3.4 | 56 | 1.000 (0.995-1.005) | 0.998 (0.987-1.011) | 1.001 (0.987-1.023) |
| Umu | 3.4-4.2 | 80 | 0.998 (0.926-1.057) | 1.002 (0.964-1.053) | 1.004 (0.988-1.033) |
| Umu | 4.2-11 | 25 | - | - | 1.001 (0.992-1.031) |
| Utau | 0-0.5 | 20 | 1.000 (1.000-1.000) | - | 1.000 (1.000-1.000) |
| Utau | 0.5-2.0 | 60 | 1.000 (0.998-1.001) | 1.000 (1.000-1.001) | 1.000 (1.000-1.015) |
| Utau | 2.0-3.4 | 56 | 1.000 (0.995-1.004) | 1.003 (0.987-1.012) | 1.002 (0.977-1.014) |
| Utau | 3.4-4.2 | 80 | 1.003 (0.938-1.089) | 0.983 (0.907-1.040) | 0.991 (0.973-1.029) |
| Utau | 4.2-11 | 25 | - | - | 0.999 (0.989-1.005) |
coarse rows vs the lxplus job rows: 312 compared, u2_min bit-identical at 7
Ue: last sensitive mass 3.66 GeV (published 3.66); best |U|^2 1.11e-08 at 0.38 GeV (published 1.11e-08 at 0.38)
Umu: last sensitive mass 3.675 GeV (published 3.665); best |U|^2 1.33e-08 at 1.575 GeV (published 1.33e-08 at 1.575)
Utau: last sensitive mass 3.855 GeV (published 3.855); best |U|^2 1.48e-07 at 3.05 GeV (published 1.48e-07 at 3.05)
Sec. IV: Umu at 1.6 GeV u2_min 1.4e-08 (published row 1.4e-08; text 1.3e-8)
Sec. IV: Ue at 1.6 GeV u2_min 1.44e-08 (published row 1.44e-08; text 1.3e-8)
Sec. IV: Ue at 0.395 GeV u2_min 1.12e-08 (published row 1.12e-08; text 1.1e-8 near 0.4)
Sec. IV: Ue min below 0.5 GeV 1.11e-08 at 0.38 (published 1.11e-08; text 1.1e-8 near 0.4 GeV)
Sec. IV: Utau best 1.48e-07 at 3.05 GeV (text ~1.5e-7)

## 5. Table (step 09) against the published table

### sensitivity.csv

rebuilt `work/hnl/table/hnl/sensitivity.csv` (723 rows) vs published `src/grendel-sensitivity/curves/data/grendel_v2_3/hnl/sensitivity.csv` (723 rows)
row order (flavour, mass) identical: True
matched rows: 723 (coarse 312, fine 411)

| column | text-identical (all) | coarse | fine | value-identical (all) |
|---|---|---|---|---|
| u2_min | 217/723 | 79/312 | 138/411 | 217/723 |
| u2_max | 302/723 | 164/312 | 138/411 | 302/723 |
| u2_min_open | 723/723 | 312/312 | 411/411 | 723/723 |
| u2_max_open | 723/723 | 312/312 | 411/411 | 723/723 |
| peak_N | 0/723 | 0/312 | 0/411 | 0/723 |
| peak_u2 | 677/723 | 304/312 | 373/411 | 677/723 |
| has_sensitivity | 721/723 | 312/312 | 409/411 | 721/723 |
| u2_min_N10 | 356/723 | 78/312 | 278/411 | 356/723 |
| u2_max_N10 | 440/723 | 162/312 | 278/411 | 440/723 |
| peak_u2_N10 | 677/723 | 304/312 | 373/411 | 677/723 |
| has_sensitivity_N10 | 723/723 | 312/312 | 411/411 | 723/723 |
| u2_min_open_N10 | 723/723 | 312/312 | 411/411 | 723/723 |
| u2_max_open_N10 | 723/723 | 312/312 | 411/411 | 723/723 |
| n_events | 723/723 | 312/312 | 411/411 | 723/723 |
| n_hits | 314/723 | 312/312 | 2/411 | 314/723 |
| n_hits_eval | 314/723 | 312/312 | 2/411 | 314/723 |
| decay_samples | 723/723 | 312/312 | 411/411 | 723/723 |
| hit_estimator | 723/723 | 312/312 | 411/411 | 723/723 |
| event_chunk | 723/723 | 312/312 | 411/411 | 723/723 |
| seed_salt | 723/723 | 312/312 | 411/411 | 723/723 |

bit-identical in every row: `u2_min_open`, `u2_max_open`, `has_sensitivity_N10`, `u2_min_open_N10`, `u2_max_open_N10`, `n_events`, `decay_samples`, `hit_estimator`, `event_chunk`, `seed_salt`

max |new/published - 1| (median in brackets) over rows where both are finite; `k one-sided` = rows where only one table has a value

| column | Ue coarse | Ue fine | Umu coarse | Umu fine | Utau coarse | Utau fine |
|---|---|---|---|---|---|---|
| u2_min | 8.06e-03 (8.2e-13) | 2.06e-01 (1.5e-03) | 2.69e-03 (1.3e-12) | 7.35e-02 (1.2e-03), 2 one-sided | 1.80e-03 (4.2e-11) | 8.92e-02 (1.9e-03) |
| u2_max | 5.83e-03 (1.8e-09) | 1.40e-01 (5.1e-03) | 2.07e-03 (2.1e-09) | 5.30e-02 (5.0e-03), 2 one-sided | 1.46e-03 (1.6e-09) | 9.28e-02 (7.8e-03) |
| peak_N | 9.58e-03 (2.1e-11) | 2.73e-02 (3.7e-03) | 7.80e-03 (2.8e-11) | 3.33e-02 (4.1e-03) | 1.48e-02 (2.0e-09) | 2.87e-02 (9.0e-03) |
| peak_u2 | 2.22e-16 (0.0e+00) | 1.20e-01 (0.0e+00) | 2.22e-16 (0.0e+00) | 1.36e-01 (0.0e+00) | 2.22e-16 (0.0e+00) | 1.36e-01 (0.0e+00) |
| u2_min_N10 | 7.53e-03 (8.8e-13) | 4.17e-02 (8.2e-04) | 2.55e-03 (1.3e-12) | 1.24e-02 (8.6e-04) | 1.43e-03 (4.4e-11) | 5.38e-02 (1.2e-03) |
| u2_max_N10 | 6.94e-03 (2.1e-09) | 2.19e-01 (5.2e-03) | 2.07e-03 (2.3e-09) | 2.15e-01 (4.5e-03) | 2.43e-03 (2.5e-09) | 2.11e-02 (4.9e-03) |
| peak_u2_N10 | 2.22e-16 (0.0e+00) | 1.20e-01 (0.0e+00) | 2.22e-16 (0.0e+00) | 1.36e-01 (0.0e+00) | 2.22e-16 (0.0e+00) | 1.36e-01 (0.0e+00) |
| n_events | 0.00e+00 (0.0e+00) | 0.00e+00 (0.0e+00) | 0.00e+00 (0.0e+00) | 0.00e+00 (0.0e+00) | 0.00e+00 (0.0e+00) | 0.00e+00 (0.0e+00) |
| n_hits | 0.00e+00 (0.0e+00) | 8.00e-03 (2.0e-03) | 0.00e+00 (0.0e+00) | 9.55e-03 (2.6e-03) | 0.00e+00 (0.0e+00) | 1.30e-02 (2.5e-03) |
| n_hits_eval | 0.00e+00 (0.0e+00) | 8.00e-03 (2.0e-03) | 0.00e+00 (0.0e+00) | 9.55e-03 (2.6e-03) | 0.00e+00 (0.0e+00) | 1.30e-02 (2.5e-03) |

verdict changes: 2

| flag | flavour | mass | block | rebuilt | published | peak N rebuilt | peak N published |
|---|---|---|---|---|---|---|---|
| has_sensitivity | Umu | 3.67 | fine | True | False | 3.073 | 2.989 |
| has_sensitivity | Umu | 3.675 | fine | True | False | 3.022 | 2.993 |

### sensitivity_nsig10.csv

rebuilt `work/hnl/table/hnl/sensitivity_nsig10.csv` (723 rows) vs published `src/grendel-sensitivity/curves/data/grendel_v2_3/hnl/sensitivity_nsig10.csv` (723 rows)
row order (flavour, mass) identical: True
matched rows: 723 (coarse 312, fine 411)

| column | text-identical (all) | coarse | fine | value-identical (all) |
|---|---|---|---|---|
| u2_min | 354/723 | 76/312 | 278/411 | 354/723 |
| u2_max | 440/723 | 162/312 | 278/411 | 440/723 |
| u2_min_open | 723/723 | 312/312 | 411/411 | 723/723 |
| u2_max_open | 723/723 | 312/312 | 411/411 | 723/723 |
| peak_N | 0/723 | 0/312 | 0/411 | 0/723 |
| peak_u2 | 447/723 | 152/312 | 295/411 | 447/723 |
| has_sensitivity | 723/723 | 312/312 | 411/411 | 723/723 |
| n_events | 723/723 | 312/312 | 411/411 | 723/723 |
| n_hits | 314/723 | 312/312 | 2/411 | 314/723 |
| n_hits_eval | 314/723 | 312/312 | 2/411 | 314/723 |
| decay_samples | 723/723 | 312/312 | 411/411 | 723/723 |
| hit_estimator | 723/723 | 312/312 | 411/411 | 723/723 |
| event_chunk | 723/723 | 312/312 | 411/411 | 723/723 |
| seed_salt | 723/723 | 312/312 | 411/411 | 723/723 |

bit-identical in every row: `u2_min_open`, `u2_max_open`, `has_sensitivity`, `n_events`, `decay_samples`, `hit_estimator`, `event_chunk`, `seed_salt`

max |new/published - 1| (median in brackets) over rows where both are finite; `k one-sided` = rows where only one table has a value

| column | Ue coarse | Ue fine | Umu coarse | Umu fine | Utau coarse | Utau fine |
|---|---|---|---|---|---|---|
| u2_min | 7.53e-03 (8.8e-13) | 4.17e-02 (8.2e-04) | 2.55e-03 (1.3e-12) | 1.24e-02 (8.6e-04) | 1.43e-03 (4.4e-11) | 5.38e-02 (1.2e-03) |
| u2_max | 6.94e-03 (2.1e-09) | 2.19e-01 (5.2e-03) | 2.07e-03 (2.3e-09) | 2.15e-01 (4.5e-03) | 2.43e-03 (2.5e-09) | 2.11e-02 (4.9e-03) |
| peak_N | 9.58e-03 (2.1e-11) | 2.73e-02 (3.7e-03) | 7.80e-03 (2.8e-11) | 3.33e-02 (4.1e-03) | 1.48e-02 (2.0e-09) | 2.87e-02 (9.0e-03) |
| peak_u2 | 6.39e-13 (0.0e+00) | 1.20e-01 (0.0e+00) | 6.76e-13 (0.0e+00) | 1.36e-01 (0.0e+00) | 6.76e-13 (3.3e-16) | 1.36e-01 (0.0e+00) |
| n_events | 0.00e+00 (0.0e+00) | 0.00e+00 (0.0e+00) | 0.00e+00 (0.0e+00) | 0.00e+00 (0.0e+00) | 0.00e+00 (0.0e+00) | 0.00e+00 (0.0e+00) |
| n_hits | 0.00e+00 (0.0e+00) | 8.00e-03 (2.0e-03) | 0.00e+00 (0.0e+00) | 9.55e-03 (2.6e-03) | 0.00e+00 (0.0e+00) | 1.30e-02 (2.5e-03) |
| n_hits_eval | 0.00e+00 (0.0e+00) | 8.00e-03 (2.0e-03) | 0.00e+00 (0.0e+00) | 9.55e-03 (2.6e-03) | 0.00e+00 (0.0e+00) | 1.30e-02 (2.5e-03) |

verdict changes: none (has_sensitivity, u2_min_open, u2_max_open)


