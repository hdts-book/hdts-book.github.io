# Production Monte Carlo driver using foreach + doParallel + %dopar% only.
# No future/parLapply backend and no resume logic.  The production output
# directory is deleted at the beginning of a run when clean_output_at_start=TRUE.

WORKER_EXPORTS_043 <- c(
  "set_rng_seed", "largest_singular_value", "spectral_radius",
  "truth_seed", "data_seed", "make_block_sizes", "make_A_true",
  "simulate_var1", "forecast_errors", "support_counts", "metrics_one",
  "assert_method_set", "fit_method_family_043", "method_replication_043"
)

SUMMARY_METRICS_043 <- c(
  "frob", "frob_sq", "msfe", "rmfe",
  "hamming", "precision", "recall", "F1", "specificity", "TPR", "FPR",
  "TP", "FP", "FN", "TN", "support_true", "support_hat",
  "off_hamming", "off_precision", "off_recall", "off_F1",
  "off_TP", "off_FP", "off_FN", "off_TN",
  "support_true_offdiag", "support_hat_offdiag",
  "rho_true", "rho_hat", "max_abs",
  "lambda_lasso", "lambda_adaptive", "scaled_lambda0",
  "base_lasso_support", "base_adaptive_support",
  "base_lasso_support_offdiag", "base_adaptive_support_offdiag",
  "threshold_factor", "runtime_fit_seconds"
)

summarize_cell_results <- function(z) {
  id_cols <- c("dgp", "d", "T_sample", "tuning", "method")
  key <- interaction(z[id_cols], drop = TRUE, lex.order = TRUE)
  pieces <- lapply(split(z, key), function(a) {
    r <- a[1, id_cols, drop = FALSE]
    r$nrep_eff <- nrow(a)
    for (v in SUMMARY_METRICS_043) {
      x <- a[[v]]
      fin <- is.finite(x)
      r[[paste0(v, "_mean")]] <- if (any(fin)) mean(x[fin]) else NA_real_
      r[[paste0(v, "_median")]] <- if (any(fin)) stats::median(x[fin]) else NA_real_
      r[[paste0(v, "_sd")]] <- if (sum(fin) >= 2L) stats::sd(x[fin]) else NA_real_
    }
    r
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

audit_truth_diagonals_043 <- function(cfg) {
  rows <- list(); kk <- 1L
  for (dgp in cfg$dgps) {
    for (idx_d in seq_along(cfg$d_grid)) {
      d <- cfg$d_grid[idx_d]
      A <- make_A_true(dgp, d, idx_d, cfg)
      rows[[kk]] <- data.frame(
        dgp = dgp, d = d,
        min_abs_diag = min(abs(diag(A))),
        n_zero_diag = sum(abs(diag(A)) <= cfg$support_tol),
        stringsAsFactors = FALSE
      )
      kk <- kk + 1L
    }
  }
  z <- do.call(rbind, rows)
  if (any(z$n_zero_diag != 0L)) stop("At least one DGP truth has a structural zero on the diagonal.")
  z
}

run_preflight_043 <- function(cfg, method_order) {
  dgp <- cfg$dgps[1]
  idx_d <- 1L
  idx_T <- 1L
  rr <- 1L
  d <- cfg$d_grid[idx_d]
  A <- make_A_true(dgp, d, idx_d, cfg)
  sd <- data_seed(cfg$master_seed, dgp, idx_d, idx_T, rr)
  z <- method_replication_043(dgp, idx_d, idx_T, rr, A, sd, cfg, method_order)
  if (nrow(z) != length(cfg$tunings) * length(method_order)) {
    stop("Preflight returned an unexpected number of method rows.")
  }

  # Critical 0.4.3 invariants.
  dat <- simulate_var1(A, cfg$T_grid[idx_T], sd, cfg)
  b <- sparseVAR::sVAR_adalasso_fista(
    Yt = dat$sample, p = cfg$p, fold = cfg$fold,
    penalize_diag = cfg$penalize_diag,
    standardize = cfg$standardize,
    whiten = cfg$whiten,
    lambda_rule = "theory", theory_c = cfg$theory_c,
    nlambda = cfg$nlambda, max_iter = cfg$max_iter, tol = cfg$tol,
    adaptive_zero = cfg$adaptive_zero, adaptive_cv = cfg$adaptive_cv
  )
  if (!identical(b$whiten, "none")) stop("Preflight: base fit is not whiten='none'.")
  if (!identical(b$penalize_diag, FALSE)) stop("Preflight: base fit penalizes the diagonal.")
  if (!is.null(b$SigmaInv)) stop("Preflight: SigmaInv must be NULL for whiten='none'.")

  Lt <- sparseVAR::sVAR_threshold(
    b, target = "lasso", factor = cfg$threshold_factor,
    penalize_diag = FALSE, p = cfg$p
  )
  At <- sparseVAR::sVAR_threshold(
    b, target = "adaptive", factor = cfg$threshold_factor,
    penalize_diag = FALSE, p = cfg$p
  )
  if (!isTRUE(all.equal(diag(Lt), diag(b$Phi_hat_lasso), tolerance = 0))) {
    stop("Preflight: L-Th changed an unpenalized diagonal coefficient.")
  }
  if (!isTRUE(all.equal(diag(At), diag(b$Phi_hat_Alasso), tolerance = 0))) {
    stop("Preflight: A-Th changed an unpenalized diagonal coefficient.")
  }
  z
}

run_methodology_simulation_043 <- function(project_dir, cfg, method_order) {
  require_replication_packages(cfg)
  assert_production_profile(cfg)
  results_dir <- safe_reset_results(project_dir, cfg$output_dir, cfg$clean_output_at_start)

  # Prevent nested BLAS/OpenMP oversubscription inside the 50 PSOCK workers.
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )

  design <- make_design_manifest(cfg)
  seeds <- make_seed_manifest(cfg)
  truths <- truth_matrix_long(cfg)
  diag_audit <- audit_truth_diagonals_043(cfg)

  truth_diag <- unique(truths[, c(
    "dgp", "idx_d", "d", "rho_true", "support_true", "support_true_offdiag"
  )])
  design <- merge(design, truth_diag, by = c("dgp", "idx_d", "d"), all.x = TRUE, sort = FALSE)
  design <- design[order(match(design$dgp, cfg$dgps), design$idx_d, design$idx_T), ]
  rownames(design) <- NULL

  write_csv(flatten_config(cfg), file.path(results_dir, "manifests", "run_config.csv"))
  write_csv(package_manifest(cfg), file.path(results_dir, "manifests", "package_versions.csv"))
  write_csv(code_manifest(project_dir), file.path(results_dir, "manifests", "code_manifest.csv"))
  write_csv(design, file.path(results_dir, "manifests", "design_manifest.csv"))
  write_csv(seeds, file.path(results_dir, "manifests", "seed_manifest.csv"))
  write_csv(truths, file.path(results_dir, "manifests", "truth_matrices.csv"))
  write_csv(diag_audit, file.path(results_dir, "manifests", "truth_diagonal_audit.csv"))

  if (isTRUE(cfg$run_preflight)) {
    cat("[preflight] package profile, diagonal convention, and first design/replication...\n")
    pf <- run_preflight_043(cfg, method_order)
    write_csv(pf, file.path(results_dir, "manifests", "preflight_result.csv"))
    cat("[preflight] PASS\n")
  }

  raw_path <- file.path(results_dir, "raw", "method_results.csv")
  if (file.exists(raw_path)) unlink(raw_path)
  cell_summaries <- vector("list", nrow(design))

  cat(sprintf(
    "[run] sparseVAR=%s | profile=no-diag/none | threshold=%.3f | cells=%d | reps/cell=%d | methods=%d | tunings=%d | workers=%d\n",
    as.character(utils::packageVersion("sparseVAR")), cfg$threshold_factor,
    nrow(design), cfg$nrep, length(method_order), length(cfg$tunings), cfg$n.cl
  ))

  cl <- parallel::makeCluster(as.integer(cfg$n.cl))
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    foreach::registerDoSEQ()
  }, add = TRUE)

  # Export helper functions exactly once.  Foreach then auto-exports only the
  # cell-specific data/arguments, avoiding the duplicate-export warnings seen
  # in the earlier diagnostic script.
  parallel::clusterExport(cl, WORKER_EXPORTS_043, envir = .GlobalEnv)
  parallel::clusterEvalQ(cl, {
    Sys.setenv(
      OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
    )
    suppressPackageStartupMessages({
      library(sparseVAR); library(MASS); library(igraph)
    })
    NULL
  })
  doParallel::registerDoParallel(cl)

  t_start <- proc.time()[3L]
  wrote_header <- FALSE
  total_rows_written <- 0L

  for (ii in seq_len(nrow(design))) {
    cc <- design[ii, , drop = FALSE]
    dgp <- as.character(cc$dgp)
    idx_d <- as.integer(cc$idx_d)
    idx_T <- as.integer(cc$idx_T)
    d <- as.integer(cc$d)
    TT <- as.integer(cc$T_sample)
    Atrue <- make_A_true(dgp, d, idx_d, cfg)
    rep_seeds <- vapply(
      seq_len(cfg$nrep),
      function(rr) data_seed(cfg$master_seed, dgp, idx_d, idx_T, rr),
      integer(1)
    )

    cat(sprintf("[cell %3d/%3d] %s d=%d T=%d ... ", ii, nrow(design), dgp, d, TT))
    flush.console()
    t0 <- proc.time()[3L]

    cell_result <- foreach::foreach(
      rr = seq_len(cfg$nrep),
      .combine = rbind,
      .inorder = TRUE,
      .errorhandling = "stop",
      .packages = character(0),
      .noexport = WORKER_EXPORTS_043
    ) %dopar% {
      method_replication_043(
        dgp = dgp, idx_d = idx_d, idx_T = idx_T, rr = rr,
        Atrue = Atrue, data_seed_value = rep_seeds[rr],
        cfg = cfg, method_order = method_order
      )
    }

    cell_result <- cell_result[
      order(cell_result$rep,
            match(cell_result$tuning, cfg$tunings),
            match(cell_result$method, method_order)),
      , drop = FALSE
    ]
    rownames(cell_result) <- NULL

    expected_rows <- cfg$nrep * length(cfg$tunings) * length(method_order)
    if (nrow(cell_result) != expected_rows) {
      stop("Unexpected row count for completed cell: ", nrow(cell_result), " != ", expected_rows)
    }
    if (length(unique(cell_result$data_seed)) != cfg$nrep) stop("Unexpected seed count in completed cell.")

    append_csv(cell_result, raw_path, write_header = !wrote_header)
    wrote_header <- TRUE
    total_rows_written <- total_rows_written + nrow(cell_result)
    cell_summaries[[ii]] <- summarize_cell_results(cell_result)

    elapsed <- proc.time()[3L] - t0
    cat(sprintf("done (%.1f min)\n", elapsed / 60))
  }

  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  on.exit(NULL, add = FALSE)

  cell_perf <- do.call(rbind, cell_summaries)
  cell_perf <- cell_perf[
    order(match(cell_perf$tuning, cfg$tunings), match(cell_perf$dgp, cfg$dgps),
          cell_perf$d, cell_perf$T_sample, match(cell_perf$method, method_order)),
    , drop = FALSE
  ]
  rownames(cell_perf) <- NULL
  write_csv(cell_perf, file.path(results_dir, "summary", "cell_performance.csv"))

  expected_raw_rows <- nrow(design) * cfg$nrep * length(cfg$tunings) * length(method_order)
  expected_cell_summary_rows <- nrow(design) * length(cfg$tunings) * length(method_order)
  if (total_rows_written != expected_raw_rows) {
    stop("Final raw-row audit failed: ", total_rows_written, " != ", expected_raw_rows)
  }
  if (nrow(cell_perf) != expected_cell_summary_rows) {
    stop("Final cell-summary audit failed: ", nrow(cell_perf), " != ", expected_cell_summary_rows)
  }

  runtime <- data.frame(
    total_elapsed_seconds = proc.time()[3L] - t_start,
    total_elapsed_hours = (proc.time()[3L] - t_start) / 3600,
    completed_design_cells = nrow(design),
    nrep = cfg$nrep,
    workers = cfg$n.cl,
    stringsAsFactors = FALSE
  )
  write_csv(runtime, file.path(results_dir, "manifests", "runtime_summary.csv"))

  integrity <- data.frame(
    item = c("raw_rows", "cell_performance_rows", "seed_manifest_rows", "design_cells"),
    expected = c(expected_raw_rows, expected_cell_summary_rows, nrow(design) * cfg$nrep, nrow(design)),
    observed = c(total_rows_written, nrow(cell_perf), nrow(seeds), nrow(design)),
    stringsAsFactors = FALSE
  )
  integrity$pass <- integrity$expected == integrity$observed
  write_csv(integrity, file.path(results_dir, "manifests", "output_integrity.csv"))
  if (!all(integrity$pass)) stop("Final output-integrity audit failed.")

  list(results_dir = results_dir, cell_performance = cell_perf)
}
