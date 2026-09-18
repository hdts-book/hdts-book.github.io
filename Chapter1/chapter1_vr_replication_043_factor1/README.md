# HDTS Chapter 1 VR replication — sparseVAR 0.4.3

This directory contains the public/book independent-design VR and VR(SNR=2) replication corresponding to the Chapter 1 sparse-estimator comparison.

## Fixed production specification

- `standardize = TRUE`
- `penalize_diag = FALSE`
- `whiten = "none"`
- `updateSigma = FALSE`
- `threshold_factor = 1.0`
- `adaptive_zero = "finite"`
- `adaptive_cv = "nested"`
- theoretical adaptive convention: `lambda1 = lambda2`
- workers: 50 PSOCK workers
- replications: 100 per `(DGP,d,T)`
- SNR settings: 1 and 2

The threshold multiplier is **not tuned**. Thresholded estimators use the canonical book rule `eta=lambda`. Unpenalized diagonal coefficients are preserved by thresholding.

## Reproducibility

The truth and data seed formulas are the same as in the VAR replication. The RNG is explicitly pinned to:

```r
RNGversion("4.3.3")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
```

with `master_seed=20260326` and `truth_seed_base=20270000`. Each `(DGP,d,T,rep)` receives its deterministic seed before data generation, so parallel scheduling does not define the random sample. SNR=1 and SNR=2 deliberately reuse the same random `X` and innovation draws for each replication; only the signal scale changes.

## Design

- DGP1--DGP3
- `d = 5,10,15,20,25,30`
- `T = 100,200,300,500,1000,2000`
- SNR = 1, 2
- theoretical and CV tuning
- 11 methods
- independent rows of `X` with covariance equal to the stationary covariance of the corresponding VAR benchmark
- innovation covariance `Sigma_e = I_d`

The raw output contains `108 * 100 * 2 * 2 * 11 = 475,200` rows.

To preserve the earlier factor-0.5 diagnostic run, this code writes to `results_vr_043_factor1/`.

## Running

Install the bundled package if needed:

```r
install.packages("sparseVAR_0.4.3_fixed.tar.gz", repos = NULL, type = "source")
```

Run:

```r
source("run_preflight_only.R")
source("run_all.R")
```

## Main outputs

- `results_vr_043_factor1/raw/vr_method_results.csv`
- `results_vr_043_factor1/summary/cell_performance.csv`
- `results_vr_043_factor1/ranks/cell_ranks.csv`
- `results_vr_043_factor1/ranks/rank_table_overall_wide.csv`
- `results_vr_043_factor1/manifests/` for seeds, configuration, package versions, code hashes, and integrity checks
