# Optimal dose regimens

For 10-day and 14-day dose regimens, we consider every dose regimen that satisfies the following criteria:

- Each daily dose comprises one or more 2.5 mg units;

- Each daily dose is at least as large as the previous dose;

- The total dose over the 10 or 14 days is 5 mg/kg; and

- We assume a body weight of 60kg, so the total dose is 300 mg (120 x 2.5 mg).

For each regimen we:

- Draw 1000 samples from the posterior (obtained by fitting the model to the ascending-dose study);

- Sample an individual random effects vector for each sample;

- Run the RBC model for each sample; and

- Record the maximum daily drop in Hb (g/dL) for each sample.

As a post-processing step, we select one or more thresholds for the maximum daily Hb drop (g/dL) and, for each dose regimen, calculate the % of individuals who exceed this threshold.

The optimal regimen(s) for a given threshold are those for which the fewest individuals exceed the threshold.

## Prerequisites

The following files must already exist:

- `../Data/RBC_model_data.RData`; and
- `../Rout/pop_fit_free_weights_cmdstan_max_delay_9_job3.rds`.

## Generate candidate dose regimens

The `candidate-dose-regimens.R` script generates all candidate dose regimens for 10-day and 14-day treatments:

```sh
./possible-dose-regimens.R
```

This script saves candidate dose regimens to the following files:

- `combinations-10-days.csv` (78,796 candidates, ~2 MB); and
- `combinations-14-days.csv` (3,648,057 candidates, ~120 MB).

## Calculate maximum daily Hb drops

**Important:** we divide the 14-day candidate regimens into 20 chunks, and process each chunk separately.
Even so, this can require a large amount of RAM, and took ~12 hours to complete on a 32-CPU virtual machine.

```sh
for _ in {1..21}; do
    ./calculate-max-daily-Hb-drops.R
done
```

This script saves results to the following 21 files:

- `evaluation-10-days-part-1.rds`;
- `evaluation-14-days-part-1.rds`;
- `evaluation-14-days-part-2.rds`;
- `evaluation-14-days-part-3.rds`;
- ...
- `evaluation-14-days-part-19.rds`; and
- `evaluation-14-days-part-20.rds`.

Note: if you intend to run this script on Debian or Ubuntu, you can use the following script to ensure that all of the necessary packages and libraries are installed:

```sh
./vm_setup.sh
```

## Identify optimal regimens for daily Hb drop thresholds

```sh
./identify-optimal-regimens.R
```

This script saves results to the following files:

- `optimal-10-day-regimens-matrix.rds`
- `optimal-10-day-regimens-table.rds`
- `optimal-14-day-regimens-matrix.rds`
- `optimal-14-day-regimens-table.rds`
- `evaluate-thresholds-results.rds`

It also creates the following plots:

- `individuals-exceeding-threshold.png`;
- `individuals-exceeding-threshold-num-optimal.png`; and
- `individuals-exceeding-threshold-optimal-regimens.png`.

## Identify near-optimal dose regimens

```sh
./evaluate-regimen-goodness.R
```

This script creates the following plots:

- `regimens-1gdl-drop-histogram.png`; and
- `regimens-1gdl-drop-intervals.png`

## Compare effective dose curves

The `compare-effective-dose-curves.R` script plots the effective dose curves for optimal 10-day and 14-day dose regimens against those of the regimens that were administered in the ascending-dose study:

```sh
compare-effective-dose-curves.R
```

This script creates the following plot:

- `compare-effective-dose-curves.png`
