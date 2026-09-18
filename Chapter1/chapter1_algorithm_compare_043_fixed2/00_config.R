# HDTS Chapter 1 optimization-backend comparison.
# Public/book profile aligned with sparseVAR 0.4.3 methodology simulations.

CFG <- list(
  project_name = "HDTS Chapter 1 solver comparison 0.4.3",
  required_sparseVAR = "0.4.3",
  strict_package_version = TRUE,

  # Design -------------------------------------------------------------------
  p = 1L,
  dgps = c("dgp1", "dgp2", "dgp3"),
  d_grid = c(5L, 10L, 15L, 20L, 25L, 30L),
  T_grid = c(100L, 200L, 300L, 500L, 1000L, 2000L),
  burn = 200L,
  T_test = 100L,

  # Monte Carlo --------------------------------------------------------------
  nrep = 100L,
  n.cl = 50L,

  # Reproducibility: identical formulas to the final VAR methodology run -----
  master_seed = 20260326L,
  truth_seed_base = 20270000L,
  rng_version = "4.3.3",
  rng_kind = "Mersenne-Twister",
  rng_normal_kind = "Inversion",
  rng_sample_kind = "Rejection",

  # DGP ----------------------------------------------------------------------
  alpha = 0.90,
  dgp1_prop = 0.10,
  dgp2_a1 = 0.50,
  dgp2_a2 = 0.90,
  dgp2_b1 = 0.20,
  dgp2_b2 = 0.40,
  dgp3_K = 2L,
  dgp3_p_in = 0.70,
  dgp3_p_out = 0.10,

  # Common estimator specification ------------------------------------------
  # All three solvers receive the SAME standardized regression, SAME penalty
  # weights, SAME lambda, and the unweighted book loss (whiten='none').
  standardize = TRUE,
  penalize_diag = FALSE,
  whiten = "none",
  theory_c = 1.0,
  adaptive_zero = "finite",
  adaptive_eps = 1e-8,
  adaptive_cap = 1e6,

  # Tuning -------------------------------------------------------------------
  # theory: exact common lambda = theory_c * sqrt(log(d^2)/N).
  # cv: common lambda1/lambda2 are selected ONCE by sparseVAR 0.4.3 FISTA
  #     using blocked/nested CV; all solvers are then fit at those same values.
  tunings = c("theory", "cv"),
  fold = 5L,
  nlambda = 80L,
  lambda_min_ratio = 1e-4,
  adaptive_cv = "nested",
  reference_max_iter = 1500L,
  reference_tol = 1e-6,

  # Final solver fits ---------------------------------------------------------
  fista_max_iter = 5000L,
  fista_tol = 1e-8,
  admm_max_iter = 5000L,
  admm_tol = 1e-7,
  admm_rho = 1.0,
  glmnet_thresh = 1e-10,
  glmnet_maxit = 1000000L,

  # Common adaptive pilot used by ALL three stage-2 solvers ------------------
  pilot_max_iter = 10000L,
  pilot_tol = 1e-10,

  # Numerical audit -----------------------------------------------------------
  support_tol = 1e-8,
  # Numerical, not statistical, zero tolerance used to canonicalize solver
  # outputs before KKT/support diagnostics.  This removes ADMM primal
  # residuals at coordinates whose proximal z-update is numerically zero.
  solver_zero_tol = 1e-5,
  preflight_max_rel_objective_spread = 5e-4,
  preflight_max_kkt = 5e-3,

  # Output -------------------------------------------------------------------
  output_dir = "results_algo_043",
  clean_output_at_start = TRUE,
  run_preflight = TRUE
)

SOLVER_ORDER <- c("glmnet", "admm", "fista")
FAMILY_ORDER <- c("lasso", "adaptive")

stopifnot(CFG$p == 1L)
stopifnot(identical(CFG$n.cl, 50L))
stopifnot(CFG$nrep == 100L)
stopifnot(identical(CFG$penalize_diag, FALSE))
stopifnot(identical(CFG$whiten, "none"))
stopifnot(isTRUE(CFG$standardize))
stopifnot(CFG$adaptive_cv == "nested")
stopifnot(CFG$solver_zero_tol > CFG$support_tol)
stopifnot(length(SOLVER_ORDER) == 3L, length(FAMILY_ORDER) == 2L)
