# Differences from `/HDTS/code/ch1.7-algo`

The final rerun changes the comparison design, not the DGP grid.

- Same DGP1--DGP3, d/T grids, burn/test lengths, and seed formula.
- Same Monte Carlo size: 100 per design cell; 50 workers.
- Adds `RNGversion("4.3.3")` for public reproducibility.
- Replaces residual covariance whitening by the book loss (`whiten="none"`).
- Keeps own-lag diagonals unpenalized consistently across all three solvers.
- Uses one common book-scale lambda under theory instead of solver-specific lambda scales.
- Uses common adaptive weights so the adaptive comparison isolates the stage-2 optimizer.
- CV is reference-tuned once and then held fixed across solvers; solver-specific CV paths are no longer compared as if they were optimizer effects.
- Adds objective/KKT/coefficient-agreement diagnostics.


- Solver outputs use a common `1e-5` numerical-zero canonicalization before KKT/support diagnostics. This addresses the ADMM backend's small primal residuals at proximal-zero coordinates and is distinct from statistical thresholding.
