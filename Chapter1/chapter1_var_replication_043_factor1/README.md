# HDTS Chapter 1 sparse-VAR replication — sparseVAR 0.4.3

This directory contains the public/book replication code for the Chapter 1 sparse-VAR estimator comparison. The DGPs, design grid, deterministic seed formulas, theory/CV comparison, Monte Carlo size, and 50-worker PSOCK/foreach backend are inherited from the preceding replication code, while the estimator specification is aligned with the Chapter 1 theory.

## Fixed production specification

```r
standardize      = TRUE
penalize_diag    = FALSE
whiten           = "none"
threshold_factor = 1.0
adaptive_zero    = "finite"
adaptive_cv      = "nested"
```

The threshold is **not tuned**. For the thresholded estimators we use the canonical book rule

```text
eta = lambda
```

so that the simulation matches the definition and support bounds for the thresholded Lasso in Chapter 1. Structurally active diagonal self-lags are not penalized and are not thresholded.

The unweighted quadratic loss (`whiten="none"`) matches the estimator analyzed in the chapter. Innovations in the VAR simulation have covariance `I_d`.

## Reproducibility

Randomness is assigned deterministically to each object and replication.

- master data-seed base: `20260326`
- truth-seed base: `20270000`
- RNG compatibility version: `RNGversion("4.3.3")`
- RNG kind: `Mersenne-Twister`
- normal kind: `Inversion`
- sample kind: `Rejection`

Each `(DGP,d,T,rep)` receives its own deterministic seed before data generation, so PSOCK worker scheduling does not determine the simulated sample. The code writes `run_config.csv`, `seed_manifest.csv`, `package_versions.csv`, and an MD5 `code_manifest.csv` with every production run.

## Production design

- DGPs: DGP1--DGP3
- `d`: 5, 10, 15, 20, 25, 30
- `T`: 100, 200, 300, 500, 1000, 2000
- replications/cell: 500
- tunings: theoretical and 5-fold blocked CV
- workers: 50 PSOCK workers
- methods: 11 Chapter 1 estimators

To preserve the earlier `threshold_factor=0.5` diagnostic run, this public factor-1 code writes to `results_043_factor1/`.

## Running

Install the bundled package if necessary:

```r
install.packages("sparseVAR_0.4.3_fixed.tar.gz", repos = NULL, type = "source")
```

Run the smoke test:

```r
source("run_preflight_only.R")
```

Then run the production simulation:

```r
source("run_all.R")
```

## Output

The main outputs are:

- `results_043_factor1/raw/method_results.csv`
- `results_043_factor1/summary/cell_performance.csv`
- `results_043_factor1/ranks/cell_ranks.csv`
- `results_043_factor1/ranks/rank_table_overall_wide.csv`
- `results_043_factor1/ranks/rank_table_offdiag_wide.csv`
- `results_043_factor1/manifests/` for seeds, configuration, versions, code hashes, and integrity checks

The four book criteria remain Frobenius loss, RMFE, normalized Hamming distance, and F1. Off-diagonal support metrics are retained as diagnostics because all DGPs contain structurally active diagonal coefficients.

## File roles

- `00_config.R`: fixed public production configuration
- `01_utils.R`: deterministic RNG/seed bookkeeping, manifests, I/O, profile checks
- `02_dgp.R`: DGP1--DGP3 and VAR(1) simulator
- `03_methods_metrics.R`: sparseVAR 0.4.3 method wrapper and metrics
- `04_simulation.R`: 50-worker production Monte Carlo driver
- `05_summaries_ranks.R`: cell summaries and ranking tables
- `run_preflight_only.R`: fast package/profile/DGP smoke test
- `run_all.R`: production entry point
- `sparseVAR_0.4.3_fixed.tar.gz`: package source used by the replication
