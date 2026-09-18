# Differences from the 0.4.2 replication

The DGP definitions, design grid, Monte Carlo seed formulas, theory/CV framework, and 50-worker execution design are preserved. The public 0.4.3 specification changes the statistical implementation as follows.

1. `penalize_diag=FALSE`: structurally active own-lag diagonal coefficients are unpenalized.
2. `whiten="none"`: the quadratic loss is unweighted, matching the Chapter 1 estimator.
3. `threshold_factor=1.0`: thresholding uses the canonical rule `eta=lambda`; no threshold multiplier is selected from simulation performance.
4. Thresholding preserves unpenalized diagonal coefficients.
5. Adaptive CV is nested, recomputing pilot fits and weights inside outer training folds.
6. The RNG semantics are pinned with `RNGversion("4.3.3")` in addition to fixed seed formulas and explicit `RNGkind`.
7. Off-diagonal support metrics are recorded as diagnostics.
8. The factor-1 run writes to `results_043_factor1/` so earlier diagnostic outputs are not overwritten.
