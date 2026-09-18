# Differences from the 0.4.2 VR replication

The DGP/SNR design, sample-size grid, truth/data seed formulas, and 100 replications per cell are preserved. The public 0.4.3 specification changes:

1. `penalize_diag=FALSE` throughout base, adaptive, thresholded, and scaled procedures.
2. `whiten="none"` and `updateSigma=FALSE`, matching the unweighted Chapter 1 loss.
3. `threshold_factor=1.0`, so thresholding is fixed at `eta=lambda` and is not tuned using simulation truth.
4. Adaptive CV remains nested.
5. RNG semantics are pinned by `RNGversion("4.3.3")` plus explicit `RNGkind` and deterministic per-replication seeds.
6. Off-diagonal support metrics are saved as diagnostics.
7. Output goes to `results_vr_043_factor1/` to preserve earlier outputs.
