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
