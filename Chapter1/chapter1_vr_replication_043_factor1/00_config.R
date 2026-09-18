# HDTS Chapter 1 VR / VR(SNR=2) replication using sparseVAR 0.4.3.
# Public/book production profile: standardized regression, no innovation
# whitening, structurally active diagonal coefficients unpenalized, and the
# canonical threshold eta=lambda (threshold_factor=1).

CFG <- list(
  project_name = "HDTS Chapter 1 VR methodology replication 0.4.3",
  required_sparseVAR = "0.4.3",
  strict_package_version = TRUE,

  # Common design grid -------------------------------------------------------
  dgps = c("dgp1", "dgp2", "dgp3"),
  d_grid = c(5L, 10L, 15L, 20L, 25L, 30L),
  T_grid = c(100L, 200L, 300L, 500L, 1000L, 2000L),
  T_test = 100L,

  # VR robustness settings --------------------------------------------------
  # SNR=1 is the main VR benchmark; SNR=2 is the stronger-signal robustness run.
  snr_grid = c(1, 2),

  # Monte Carlo -------------------------------------------------------------
  # Preserved from the completed 0.4.2 VR replication.
  nrep = 100L,
  n.cl = 50L,
  master_seed = 20260326L,
  truth_seed_base = 20270000L,
  rng_version = "4.3.3",
  rng_kind = "Mersenne-Twister",
  rng_normal_kind = "Inversion",
  rng_sample_kind = "Rejection",

  # DGP ---------------------------------------------------------------------
  alpha = 0.90,
  dgp1_prop = 0.10,
  dgp2_a1 = 0.50,
  dgp2_a2 = 0.90,
  dgp2_b1 = 0.20,
  dgp2_b2 = 0.40,
  dgp3_K = 2L,
  dgp3_p_in = 0.70,
  dgp3_p_out = 0.10,

  # sparseVAR 0.4.3 / book-theory profile ----------------------------------
  # The Chapter-1 theory does not use Sigma^{-1}-weighted loss.  Hence the
  # production benchmark uses the naive/unweighted quadratic loss.
  standardize = TRUE,
  penalize_diag = FALSE,
  whiten = "none",               # one of: none, diag, full
  updateSigma = FALSE,            # must be FALSE when whiten == "none"

  # Tuning ------------------------------------------------------------------
  tunings = c("theory", "cv"),
  fold = 5L,
  nlambda = 80L,
  lambda_min_ratio = 1e-4,
  theory_c = 1.0,
  max_iter = 1500L,
  tol = 1e-6,

  # Adaptive Lasso ----------------------------------------------------------
  adaptive_zero = "finite",
  adaptive_cv = "nested",
  adaptive_eps = 1e-8,
  adaptive_cap = 1e6,

  # Theory convention retained for the book:
  # lambda1 = lambda2 = sqrt(log(q)/N), q = d^2.
  theory_same_lambda = TRUE,

  # Post-processing ---------------------------------------------------------
  # Fixed a priori to match the Chapter 1 thresholded-Lasso definition.
  threshold_factor = 1.0,
  relax_gamma = 0.50,
  support_tol = 1e-8,
  scaled_lambda0 = NULL,
  scaled_max_outer = 100L,
  scaled_tol_outer = 1e-6,

  # Output ------------------------------------------------------------------
  output_dir = "results_vr_043_factor1",
  clean_output_at_start = TRUE,
  run_preflight = TRUE
)

METHOD_ORDER <- c(
  "L", "L-Th", "L-Ref", "L-Rlx", "L-Des",
  "A", "A-Th", "A-Ref", "A-Rlx", "A-Des",
  "ScL"
)

METHOD_LABELS <- c(
  L = "Lasso",
  `L-Th` = "Lasso+Thr",
  `L-Ref` = "Lasso+Refit",
  `L-Rlx` = "Lasso+Relax",
  `L-Des` = "Lasso+DeSP",
  A = "Ada",
  `A-Th` = "Ada+Thr",
  `A-Ref` = "Ada+Refit",
  `A-Rlx` = "Ada+Relax",
  `A-Des` = "Ada+DeSP",
  ScL = "Scaled"
)

stopifnot(CFG$n.cl == 50L)
stopifnot(CFG$nrep == 100L)
stopifnot(all(CFG$snr_grid > 0))
stopifnot(CFG$adaptive_cv == "nested")
stopifnot(CFG$adaptive_zero == "finite")
stopifnot(CFG$whiten %in% c("none", "diag", "full"))
stopifnot(CFG$penalize_diag == FALSE)
stopifnot(identical(CFG$threshold_factor, 1.0))
stopifnot(identical(CFG$rng_version, "4.3.3"))
if (CFG$whiten == "none" && isTRUE(CFG$updateSigma)) {
  stop("Inconsistent configuration: whiten='none' requires updateSigma=FALSE.")
}
if (CFG$whiten != "none" && !isTRUE(CFG$updateSigma)) {
  stop("Inconsistent configuration: diag/full whitening requires updateSigma=TRUE.")
}
