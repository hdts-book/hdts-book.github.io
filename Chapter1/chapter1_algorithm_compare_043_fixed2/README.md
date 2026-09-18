# HDTS Chapter 1: optimization-backend comparison (sparseVAR 0.4.3)

This directory reruns the Chapter 1 comparison of `glmnet`, ADMM, and FISTA under the **final book estimator specification**.

## Statistical profile

- `standardize = TRUE`
- `penalize_diag = FALSE`
- `whiten = "none"`
- innovation covariance in the DGP is exactly `I_d`
- fixed `RNGversion("4.3.3")`
- `master_seed = 20260326`, `truth_seed_base = 20270000`
- 100 replications per `(DGP,d,T)` cell
- 50 PSOCK workers

The DGPs, dimensions, sample-size grid, truth seeds, and data-seed formula are identical to the final Chapter 1 VAR methodology simulation. Thus the first 100 replications of each cell use the same simulated samples.

## What is different from the old solver comparison

The earlier `ch1.7-algo` experiment mixed solver-specific preprocessing and penalty scales. In particular it used residual covariance updates/whitening and did not put all three solvers on the same numerical objective.

This version deliberately isolates the optimizer:

1. All solvers receive the same centered/standardized `X,Y` matrices.
2. All use the book loss `(1/N)||Y-XB||_F^2` with no whitening.
3. Own-lag diagonal coefficients are unpenalized for every solver.
4. Under theoretical tuning, all solvers receive the exact same `lambda = sqrt(log(d^2)/N)`.
5. The adaptive comparison uses one common high-accuracy FISTA pilot to construct a single weight matrix; glmnet, ADMM and FISTA then solve the **same weighted-Lasso problem**.
6. The CV experiment is secondary: sparseVAR 0.4.3 FISTA with blocked/nested CV selects common `lambda1,lambda2` once, after which all three final solvers are evaluated at the same lambda values and common adaptive weights. This prevents solver-specific lambda paths from contaminating the optimizer comparison.

Because the statistical objective is now identical, large differences in Frobenius/RMFE/support metrics would indicate numerical/scaling problems rather than meaningful estimator differences. For this reason the output also reports objective gaps, KKT violations, coefficient differences to FISTA, and ADMM convergence.

## Running

Install the bundled sparseVAR source if needed:

```r
install.packages("sparseVAR_0.4.3_fixed.tar.gz", repos = NULL, type = "source")
```

Then run the fast smoke test:

```r
source("run_preflight_only.R")
```

The preflight checks that the three solvers agree on the same small theory-tuned problem. If it passes, the main book comparison can be regenerated with:

```r
source("run_theory_only.R")
```

This is the primary evidence for the optimization-backend subsection. The supplementary common-tuning CV diagnostic is:

```r
source("run_cv_only.R")
```

To run both theory and CV into one result directory:

```r
source("run_all.R")
```

## Key outputs

- `raw/solver_results.csv`
- `summary/cell_performance.csv`
- `ranks/rank_table_theory.csv`
- `ranks/rank_table_cv.csv`
- `summary/winner_map_data.csv`
- `summary/numerical_agreement.csv`
- `manifests/seed_manifest.csv`
- `manifests/run_config.csv`
- `manifests/package_versions.csv`
- `manifests/code_manifest.csv`
- `manifests/output_integrity.csv`

For a theory-only run, the expected raw row count is
`108 cells x 100 reps x 2 families x 3 solvers = 64,800`.
For a combined theory+CV run it is `129,600`.

### Rank ties

Since the three backends are now solving the same convex problem, their statistical metrics can differ only at numerical precision when all solvers converge. Cell ranks and winner labels therefore use an 8-significant-digit tolerance for tie handling. The unrounded objective, KKT, and coefficient-discrepancy diagnostics remain available in the raw output.

## Numerical-zero convention

The ADMM backend returns its quadratic-update variable rather than the proximal splitting variable. At convergence these agree up to the primal residual, but entries that are structurally zero can therefore remain at roughly `1e-7`--a few `1e-6`. For solver comparison only, all three standardized coefficient matrices are canonicalized by setting `|B_ij| <= 1e-5` to exactly zero before KKT and support diagnostics. This is a numerical cleanup, not the statistical thresholded-Lasso operation. The number of cleaned entries and the largest cleaned magnitude are saved in the raw output.

## PSOCK worker export note

All helper functions used inside the parallel solver loop are explicitly exported to PSOCK workers. In particular, `canonicalize_solver_B()` is exported because `fit_solver_B()` calls it after each solver fit. This avoids master-process-only preflight success followed by worker lookup failures in the production loop.
