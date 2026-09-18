# Production VR / VR(SNR=2) simulation for sparseVAR 0.4.3.
# The 0.4.2 SNR design, truth/data seeds, and R=100 are preserved.

`%dopar%` <- foreach::`%dopar%`

make_manifest_vr <- function(cfg) {
  z <- expand.grid(
    dgp = cfg$dgps,
    idx_d = seq_along(cfg$d_grid),
    idx_T = seq_along(cfg$T_grid),
    rep = seq_len(cfg$nrep),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  z$d <- cfg$d_grid[z$idx_d]
  z$T_sample <- cfg$T_grid[z$idx_T]
  z$truth_seed <- mapply(
    truth_seed, z$dgp, z$idx_d,
    MoreArgs = list(base = cfg$truth_seed_base)
  )
  z$data_seed <- mapply(
    data_seed, z$dgp, z$idx_d, z$idx_T, z$rep,
    MoreArgs = list(master_seed = cfg$master_seed)
  )
  z <- z[order(match(z$dgp, cfg$dgps), z$idx_d, z$idx_T, z$rep), ]
  rownames(z) <- NULL
  if (anyDuplicated(z$data_seed)) stop("Duplicate data seeds detected.")
  z
}

one_vr_task <- function(task, cfg, method_order) {
  dgp <- as.character(task$dgp)
  idx_d <- as.integer(task$idx_d)
  idx_T <- as.integer(task$idx_T)
  d <- as.integer(task$d)
  TT <- as.integer(task$T_sample)
  rr <- as.integer(task$rep)
  sd <- as.integer(task$data_seed)

  A_base <- make_A_true(dgp, d, idx_d, cfg)

  out <- vector("list", length(cfg$snr_grid) * length(cfg$tunings) * length(method_order))
  kk <- 1L

  for (snr in cfg$snr_grid) {
    # The same X/E random numbers are intentionally paired across SNR settings,
    # exactly as in the 0.4.2 VR benchmark; only the signal scale changes.
    dat <- simulate_vr(
      A_base = A_base,
      T_sample = TT,
      T_test = cfg$T_test,
      snr_target = snr,
      seed = sd,
      cfg = cfg
    )

    for (tuning in cfg$tunings) {
      fit <- fit_method_family_vr(dat$X_sample, dat$Y_sample, tuning, cfg)
      if (!setequal(names(fit$methods), method_order)) stop("Method set mismatch.")

      for (mm in method_order) {
        m <- metrics_vr(
          fit$methods[[mm]], dat$A_true,
          dat$X_test, dat$Y_test, cfg
        )
        m$rep <- rr
        m$dgp <- dgp
        m$idx_d <- idx_d
        m$d <- d
        m$idx_T <- idx_T
        m$T_sample <- TT
        m$snr <- snr
        m$tuning <- tuning
        m$method <- mm
        m$data_seed <- sd
        m$truth_seed <- as.integer(task$truth_seed)
        m$signal_scale <- dat$signal_scale
        m$lambda_lasso <- fit$base$lambda_lasso
        m$lambda_adaptive <- if (startsWith(mm, "A")) fit$base$lambda_adalasso else NA_real_
        m$scaled_lambda0 <- if (mm == "ScL") fit$scaled$lambda0 else NA_real_
        m$base_lasso_support <- sum(abs(fit$base$Phi_hat_lasso) > cfg$support_tol)
        m$base_adaptive_support <- sum(abs(fit$base$Phi_hat_Alasso) > cfg$support_tol)
        m$runtime_fit_seconds <- fit$runtime

        out[[kk]] <- m[, c(
          "rep", "dgp", "idx_d", "d", "idx_T", "T_sample",
          "snr", "tuning", "method",
          "data_seed", "truth_seed", "signal_scale",
          "frob", "rmfe", "hamming", "F1",
          "precision", "recall", "TP", "FP", "FN", "TN",
          "support_true", "support_hat",
          "off_hamming", "off_F1", "off_precision", "off_recall",
          "off_TP", "off_FP", "off_FN", "off_TN",
          "off_support_true", "off_support_hat",
          "lambda_lasso", "lambda_adaptive", "scaled_lambda0",
          "base_lasso_support", "base_adaptive_support",
          "runtime_fit_seconds"
        )]
        kk <- kk + 1L
      }
    }
  }

  ans <- do.call(rbind, out)
  core <- c("frob", "rmfe", "hamming", "F1", "off_hamming")
  if (any(!is.finite(as.matrix(ans[, core, drop = FALSE])))) {
    stop("Non-finite core metric: ", dgp, " d=", d, " T=", TT, " rep=", rr)
  }
  ans
}

run_preflight_vr_043 <- function(cfg, method_order) {
  task <- data.frame(
    dgp = cfg$dgps[1],
    idx_d = 1L,
    idx_T = 1L,
    rep = 1L,
    d = cfg$d_grid[1],
    T_sample = cfg$T_grid[1],
    truth_seed = truth_seed(cfg$dgps[1], 1L, cfg$truth_seed_base),
    data_seed = data_seed(cfg$master_seed, cfg$dgps[1], 1L, 1L, 1L),
    stringsAsFactors = FALSE
  )
  z <- one_vr_task(task, cfg, method_order)
  expected <- length(cfg$snr_grid) * length(cfg$tunings) * length(method_order)
  if (nrow(z) != expected) stop("Preflight row-count mismatch.")
  if (!all(z$off_FP >= 0 & z$off_FN >= 0)) stop("Invalid off-diagonal support metrics.")
  z
}

run_vr_simulation <- function(cfg, project_dir, method_order) {
  require_vr_packages(cfg)
  if (cfg$n.cl != 50L) stop("Production VR run requires exactly 50 PSOCK workers.")
  if (!identical(cfg$threshold_factor, 1.0)) stop("Production VR profile requires threshold_factor=1.0 (eta=lambda).")
  if (!identical(cfg$rng_version, "4.3.3")) stop("Production VR profile requires rng_version=\"4.3.3\".")

  out_root <- file.path(project_dir, cfg$output_dir)
  if (isTRUE(cfg$clean_output_at_start) && dir.exists(out_root)) {
    unlink(out_root, recursive = TRUE, force = TRUE)
  }
  for (dd in c("raw", "summary", "ranks", "manifests")) {
    dir.create(file.path(out_root, dd), recursive = TRUE, showWarnings = FALSE)
  }

  # Avoid BLAS/OpenMP oversubscription inside 50 R workers.
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )

  manifest <- make_manifest_vr(cfg)
  utils::write.csv(manifest, file.path(out_root, "manifests", "seed_manifest.csv"), row.names = FALSE)
  utils::write.csv(flatten_config(cfg), file.path(out_root, "manifests", "run_config.csv"), row.names = FALSE)
  utils::write.csv(package_manifest(), file.path(out_root, "manifests", "package_versions.csv"), row.names = FALSE)
  utils::write.csv(code_manifest(project_dir), file.path(out_root, "manifests", "code_manifest.csv"), row.names = FALSE)

  # Truth matrices, one per (DGP,d), on the base scale.
  truth_rows <- list()
  kk <- 1L
  for (g in cfg$dgps) {
    for (idx_d in seq_along(cfg$d_grid)) {
      d <- cfg$d_grid[idx_d]
      A <- make_A_true(g, d, idx_d, cfg)
      ij <- expand.grid(row = seq_len(d), col = seq_len(d))
      truth_rows[[kk]] <- data.frame(
        dgp = g,
        idx_d = idx_d,
        d = d,
        truth_seed = truth_seed(g, idx_d, cfg$truth_seed_base),
        row = ij$row,
        col = ij$col,
        A_base = A[cbind(ij$row, ij$col)],
        stringsAsFactors = FALSE
      )
      kk <- kk + 1L
    }
  }
  utils::write.csv(
    do.call(rbind, truth_rows),
    file.path(out_root, "manifests", "truth_matrices_base.csv"),
    row.names = FALSE
  )

  if (isTRUE(cfg$run_preflight)) {
    cat("[preflight] VR dgp1/d=5/T=100, SNR=1/2, theory/CV ... ")
    pf <- run_preflight_vr_043(cfg, method_order)
    utils::write.csv(pf, file.path(out_root, "manifests", "preflight_result.csv"), row.names = FALSE)
    cat("PASS\n")
  }

  cl <- parallel::makePSOCKcluster(as.integer(cfg$n.cl))
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    foreach::registerDoSEQ()
  }, add = TRUE)

  worker_exports <- c(
    "CFG", "METHOD_ORDER",
    "set_rng_seed", "largest_singular_value", "spectral_radius",
    "truth_seed", "data_seed", "make_block_sizes", "make_A_true",
    "stationary_cov_var1", "signal_scale_for_snr", "simulate_vr",
    "block_foldid", "safe_invsqrt", "prepare_vr", "fista_backend",
    "weighted_sse", "make_penalty_weights", "kkt_lambda_max", "lambda_grid",
    "cv_lambda", "adaptive_weights", "nested_adaptive_lambda",
    "backtransform", "fit_vr_fista", "vr_refit", "vr_debias",
    "vr_threshold", "vr_scaled_lasso", "fit_method_family_vr",
    "forecast_rmfe_vr", "binary_support_metrics", "metrics_vr",
    "one_vr_task"
  )
  parallel::clusterExport(cl, worker_exports, envir = .GlobalEnv)

  tasks <- split(manifest, seq_len(nrow(manifest)))
  t0 <- proc.time()[3L]

  cat(sprintf(
    "[run] sparseVAR=%s | tasks=%d | reps/cell=%d | SNR=%s | tunings=%d | methods=%d | workers=%d\n",
    as.character(utils::packageVersion("sparseVAR")),
    length(tasks), cfg$nrep, paste(cfg$snr_grid, collapse = ","),
    length(cfg$tunings), length(method_order), cfg$n.cl
  ))

  raw <- tryCatch(
    {
      foreach::foreach(
        task = tasks,
        .combine = rbind,
        .multicombine = TRUE,
        .maxcombine = 100L,
        .inorder = FALSE,
        .errorhandling = "stop",
        .packages = c("sparseVAR", "MASS", "igraph"),
        .noexport = worker_exports
      ) %dopar% {
        one_vr_task(task, CFG, METHOD_ORDER)
      }
    },
    finally = {
      try(parallel::stopCluster(cl), silent = TRUE)
      foreach::registerDoSEQ()
    }
  )
  on.exit(NULL, add = FALSE)

  elapsed <- proc.time()[3L] - t0

  raw <- raw[order(
    match(raw$dgp, cfg$dgps),
    raw$d, raw$T_sample, raw$rep,
    match(raw$snr, cfg$snr_grid),
    match(raw$tuning, cfg$tunings),
    match(raw$method, method_order)
  ), , drop = FALSE]
  rownames(raw) <- NULL

  expected <- length(cfg$dgps) * length(cfg$d_grid) * length(cfg$T_grid) *
    cfg$nrep * length(cfg$snr_grid) * length(cfg$tunings) * length(method_order)
  if (nrow(raw) != expected) stop("Raw row count mismatch: ", nrow(raw), " vs ", expected)
  if (length(unique(raw$data_seed)) != nrow(manifest)) {
    stop("Unexpected unique data-seed count in raw output.")
  }

  utils::write.csv(raw, file.path(out_root, "raw", "vr_method_results.csv"), row.names = FALSE)

  integrity <- data.frame(
    item = c(
      "raw_rows", "base_design_cells", "manifest_rows",
      "snr_settings", "tunings", "methods", "workers"
    ),
    expected = c(
      expected, 108L, 108L * cfg$nrep,
      length(cfg$snr_grid), length(cfg$tunings), length(method_order), 50L
    ),
    observed = c(
      nrow(raw),
      length(cfg$dgps) * length(cfg$d_grid) * length(cfg$T_grid),
      nrow(manifest),
      length(unique(raw$snr)), length(unique(raw$tuning)),
      length(unique(raw$method)), cfg$n.cl
    ),
    stringsAsFactors = FALSE
  )
  integrity$pass <- integrity$expected == integrity$observed
  utils::write.csv(integrity, file.path(out_root, "manifests", "output_integrity.csv"), row.names = FALSE)
  if (!all(integrity$pass)) stop("Output-integrity audit failed.")

  runtime <- data.frame(
    elapsed_seconds = elapsed,
    elapsed_hours = elapsed / 3600,
    workers = cfg$n.cl,
    nrep = cfg$nrep,
    stringsAsFactors = FALSE
  )
  utils::write.csv(runtime, file.path(out_root, "manifests", "runtime_summary.csv"), row.names = FALSE)

  raw
}
