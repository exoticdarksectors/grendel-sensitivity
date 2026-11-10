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
