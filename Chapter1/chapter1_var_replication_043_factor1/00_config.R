# HDTS Chapter 1 sparse-VAR methodology replication
# Production profile: sparseVAR 0.4.3, book-theory-aligned naive loss.
#
# This file intentionally keeps the 0.4.2 design/seed/parallel structure.
# Changing any value below defines a different experiment.

CFG <- list(
  project_name = "HDTS Chapter 1 sparse-VAR methodology replication 0.4.3",
  required_sparseVAR = "0.4.3",
  strict_package_version = TRUE,

  # Model / design grid -------------------------------------------------------
  p = 1L,
  dgps = c("dgp1", "dgp2", "dgp3"),
  d_grid = c(5L, 10L, 15L, 20L, 25L, 30L),
  T_grid = c(100L, 200L, 300L, 500L, 1000L, 2000L),
  burn = 200L,
  T_test = 100L,

  # Monte Carlo ---------------------------------------------------------------
  nrep = 500L,
  n.cl = 50L,

  # Seeds: unchanged from the 0.4.2 production replication -------------------
  master_seed = 20260326L,
  truth_seed_base = 20270000L,
  rng_version = "4.3.3",
  rng_kind = "Mersenne-Twister",
  rng_normal_kind = "Inversion",
  rng_sample_kind = "Rejection",

  # DGP parameters: unchanged -------------------------------------------------
  alpha = 0.90,
  dgp1_prop = 0.10,
  dgp2_a1 = 0.50,
  dgp2_a2 = 0.90,
  dgp2_b1 = 0.20,
  dgp2_b2 = 0.40,
  dgp3_K = 2L,
  dgp3_p_in = 0.70,
  dgp3_p_out = 0.10,

  # sparseVAR 0.4.3 public/book profile ---------------------------------------
  # Diagonal self-lags are structurally active in DGP1--DGP3, so they are
  # treated as unpenalized.  The book theory uses the unweighted RSS loss, so
  # no estimated innovation-covariance whitening is used.
  penalize_diag = FALSE,
  standardize = TRUE,
  whiten = "none",

  # Tuning --------------------------------------------------------------------
  tunings = c("theory", "cv"),
  fold = 5L,
  nlambda = 80L,
  theory_c = 1.0,
  max_iter = 1500L,
  tol = 1e-6,

  # Adaptive Lasso ------------------------------------------------------------
  adaptive_zero = "finite",
  adaptive_cv = "nested",

  # Post-processing -----------------------------------------------------------
  # Canonical book specification: threshold at the same regularization scale,
  # eta = lambda.  The threshold constant is fixed a priori and is not tuned.
  threshold_factor = 1.0,
  relax_gamma = 0.50,
  scaled_lambda0 = NULL,  # NULL -> Sun--Zhang sqrt(2 log(d*p)/N) in sparseVAR
  support_tol = 1e-8,

  # Production output ---------------------------------------------------------
  output_dir = "results_043_factor1",
  clean_output_at_start = TRUE,
  run_preflight = TRUE
)

METHOD_ORDER <- c(
  "L", "L-Th", "L-Ref", "L-Rlx", "L-Des",
  "A", "A-Th", "A-Ref", "A-Rlx", "A-Des", "ScL"
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

stopifnot(CFG$p == 1L)
stopifnot(identical(CFG$n.cl, 50L))
stopifnot(CFG$nrep == 500L)
stopifnot(identical(CFG$penalize_diag, FALSE))
stopifnot(identical(CFG$whiten, "none"))
stopifnot(isTRUE(CFG$standardize))
stopifnot(CFG$adaptive_zero %in% c("finite", "active"))
stopifnot(CFG$adaptive_cv %in% c("nested", "fixed_weights"))
stopifnot(CFG$fold >= 2L, CFG$nlambda >= 2L)
stopifnot(length(METHOD_ORDER) == 11L)
stopifnot(identical(CFG$threshold_factor, 1.0))
stopifnot(identical(CFG$rng_version, "4.3.3"))
