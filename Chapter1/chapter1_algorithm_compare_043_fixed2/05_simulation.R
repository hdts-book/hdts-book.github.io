# 50-worker deterministic Monte Carlo driver.

`%dopar%` <- foreach::`%dopar%`

WORKER_EXPORTS_ALGO <- c(
  "set_rng_seed", "largest_singular_value", "spectral_radius", "truth_seed", "data_seed",
  "make_block_sizes", "make_A_true", "simulate_var1",
  "prepare_common_regression", "make_penalty_weights", "backtransform_B", "book_theory_lambda",
  "fit_fista_B", "fit_admm_B", "fit_glmnet_B", "canonicalize_solver_B", "fit_solver_B", "adaptive_weights_common",
  "reference_tuning", "objective_book", "kkt_violation",
  "forecast_errors", "support_counts", "metrics_one", "one_algorithm_replication"
)

SUMMARY_METRICS_ALGO <- c(
  "frob", "frob_sq", "msfe", "rmfe", "hamming", "precision", "recall", "F1",
  "TP", "FP", "FN", "TN", "support_true", "support_hat",
  "off_hamming", "off_precision", "off_recall", "off_F1",
  "off_TP", "off_FP", "off_FN", "off_TN", "support_true_offdiag", "support_hat_offdiag",
  "rho_true", "rho_hat", "runtime_seconds", "objective", "objective_gap", "kkt_max",
  "n_solver_zeroed", "max_solver_zeroed_abs",
  "B_diff_to_fista", "Phi_diff_to_fista", "reference_tuning_seconds", "reference_pilot_seconds"
)

summarize_algorithm_cells <- function(raw) {
  idc <- c("tuning", "family", "solver", "dgp", "d", "T_sample")
  key <- interaction(raw[idc], drop = TRUE, lex.order = TRUE)
  out <- do.call(rbind, lapply(split(raw, key), function(a) {
    r <- a[1, idc, drop = FALSE]
    r$nrep_eff <- nrow(a)
    for (v in SUMMARY_METRICS_ALGO) {
      x <- a[[v]]; fin <- is.finite(x)
      r[[paste0(v, "_mean")]] <- if (any(fin)) mean(x[fin]) else NA_real_
      r[[paste0(v, "_median")]] <- if (any(fin)) stats::median(x[fin]) else NA_real_
      r[[paste0(v, "_sd")]] <- if (sum(fin) >= 2L) stats::sd(x[fin]) else NA_real_
    }
    r$admm_convergence_rate <- if (r$solver == "admm") mean(a$admm_converged, na.rm = TRUE) else NA_real_
    r
  }))
  rownames(out) <- NULL
  out
}

run_preflight_algorithm <- function(cfg, solver_order, family_order) {
  dgp <- "dgp1"; idx_d <- 1L; idx_T <- 1L; rr <- 1L
  A <- make_A_true(dgp, cfg$d_grid[idx_d], idx_d, cfg)
  sd <- data_seed(cfg$master_seed, dgp, idx_d, idx_T, rr)
  cfg0 <- cfg; cfg0$tunings <- "theory"
  z <- one_algorithm_replication(dgp, idx_d, idx_T, rr, A, sd, cfg0, solver_order, family_order)
  diag_cols <- c("family", "solver", "lambda_used", "objective", "objective_gap",
                 "kkt_max", "B_diff_to_fista", "n_solver_zeroed",
                 "max_solver_zeroed_abs", "admm_converged")
  print(z[, diag_cols, drop = FALSE], row.names = FALSE)
  for (fam in family_order) {
    a <- z[z$family == fam, , drop = FALSE]
    rel <- (max(a$objective) - min(a$objective)) / max(1, abs(min(a$objective)))
    if (rel > cfg$preflight_max_rel_objective_spread) {
      stop("Preflight objective mismatch in ", fam, ": relative spread=", signif(rel, 4))
    }
    if (max(a$kkt_max) > cfg$preflight_max_kkt) {
      ibad <- which.max(a$kkt_max)
      stop(
        "Preflight KKT mismatch in ", fam,
        " | solver=", a$solver[ibad],
        " | max=", signif(a$kkt_max[ibad], 4),
        " | objective_gap=", signif(a$objective_gap[ibad], 4),
        " | max_solver_zeroed_abs=", signif(a$max_solver_zeroed_abs[ibad], 4)
      )
    }
  }
  z
}

run_algorithm_simulation <- function(project_dir, cfg, solver_order, family_order) {
  require_algorithm_packages(cfg)
  results_dir <- safe_reset_results(project_dir, cfg$output_dir, cfg$clean_output_at_start)
  Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

  design <- make_design_manifest(cfg)
  seeds <- make_seed_manifest(cfg)
  write_csv(flatten_config(cfg), file.path(results_dir, "manifests", "run_config.csv"))
  write_csv(package_manifest(), file.path(results_dir, "manifests", "package_versions.csv"))
  write_csv(code_manifest(project_dir), file.path(results_dir, "manifests", "code_manifest.csv"))
  write_csv(design, file.path(results_dir, "manifests", "design_manifest.csv"))
  write_csv(seeds, file.path(results_dir, "manifests", "seed_manifest.csv"))

  if (isTRUE(cfg$run_preflight)) {
    cat("[preflight] common-objective agreement across glmnet/ADMM/FISTA ... ")
    pf <- run_preflight_algorithm(cfg, solver_order, family_order)
    write_csv(pf, file.path(results_dir, "manifests", "preflight_result.csv"))
    cat("PASS\n")
  }

  cl <- parallel::makePSOCKcluster(as.integer(cfg$n.cl))
  on.exit({ try(parallel::stopCluster(cl), silent=TRUE); foreach::registerDoSEQ() }, add=TRUE)
  parallel::clusterExport(cl, WORKER_EXPORTS_ALGO, envir=.GlobalEnv)
  parallel::clusterEvalQ(cl, {
    Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    suppressPackageStartupMessages({ library(sparseVAR); library(glmnet); library(MASS); library(igraph) })
    NULL
  })
  doParallel::registerDoParallel(cl)

  raw_path <- file.path(results_dir, "raw", "solver_results.csv")
  if (file.exists(raw_path)) unlink(raw_path)
  wrote_header <- FALSE
  summaries <- vector("list", nrow(design))
  total_rows <- 0L
  t_start <- proc.time()[3L]

  cat(sprintf("[run] cells=%d reps/cell=%d tunings=%s families=%d solvers=%d workers=%d\n",
              nrow(design), cfg$nrep, paste(cfg$tunings, collapse=","),
              length(family_order), length(solver_order), cfg$n.cl))

  for (ii in seq_len(nrow(design))) {
    cc <- design[ii, , drop=FALSE]
    dgp <- as.character(cc$dgp); idx_d <- as.integer(cc$idx_d); idx_T <- as.integer(cc$idx_T)
    d <- as.integer(cc$d); TT <- as.integer(cc$T_sample)
    A <- make_A_true(dgp, d, idx_d, cfg)
    rep_seeds <- vapply(seq_len(cfg$nrep), function(rr) data_seed(cfg$master_seed, dgp, idx_d, idx_T, rr), integer(1))
    cat(sprintf("[cell %3d/%3d] %s d=%d T=%d ... ", ii, nrow(design), dgp, d, TT)); flush.console()
    t0 <- proc.time()[3L]

    cell <- foreach::foreach(
      rr = seq_len(cfg$nrep), .combine=rbind, .inorder=TRUE,
      .errorhandling="stop", .packages=character(0), .noexport=WORKER_EXPORTS_ALGO
    ) %dopar% {
      one_algorithm_replication(dgp, idx_d, idx_T, rr, A, rep_seeds[rr], cfg, solver_order, family_order)
    }

    cell <- cell[order(cell$rep, match(cell$tuning,cfg$tunings),
                       match(cell$family,family_order), match(cell$solver,solver_order)), , drop=FALSE]
    expected <- cfg$nrep * length(cfg$tunings) * length(family_order) * length(solver_order)
    if (nrow(cell) != expected) stop("Cell row count mismatch.")
    if (length(unique(cell$data_seed)) != cfg$nrep) stop("Cell seed count mismatch.")
    append_csv(cell, raw_path, write_header=!wrote_header); wrote_header <- TRUE
    summaries[[ii]] <- summarize_algorithm_cells(cell)
    total_rows <- total_rows + nrow(cell)
    cat(sprintf("done (%.1f min)\n", (proc.time()[3L]-t0)/60))
  }

  parallel::stopCluster(cl); foreach::registerDoSEQ(); on.exit(NULL, add=FALSE)
  cell_perf <- do.call(rbind, summaries)
  cell_perf <- cell_perf[order(match(cell_perf$tuning,cfg$tunings), match(cell_perf$family,family_order),
                               match(cell_perf$dgp,cfg$dgps), cell_perf$d, cell_perf$T_sample,
                               match(cell_perf$solver,solver_order)), , drop=FALSE]
  rownames(cell_perf) <- NULL
  write_csv(cell_perf, file.path(results_dir, "summary", "cell_performance.csv"))

  expected_raw <- nrow(design)*cfg$nrep*length(cfg$tunings)*length(family_order)*length(solver_order)
  expected_cell <- nrow(design)*length(cfg$tunings)*length(family_order)*length(solver_order)
  integrity <- data.frame(
    item=c("raw_rows","cell_performance_rows","seed_manifest_rows","design_cells","workers"),
    expected=c(expected_raw,expected_cell,nrow(design)*cfg$nrep,nrow(design),cfg$n.cl),
    observed=c(total_rows,nrow(cell_perf),nrow(seeds),nrow(design),cfg$n.cl),
    stringsAsFactors=FALSE)
  integrity$pass <- integrity$expected == integrity$observed
  write_csv(integrity, file.path(results_dir,"manifests","output_integrity.csv"))
  if (!all(integrity$pass)) stop("Output integrity audit failed.")
  write_csv(data.frame(elapsed_seconds=proc.time()[3L]-t_start,
                       elapsed_hours=(proc.time()[3L]-t_start)/3600,
                       workers=cfg$n.cl,nrep=cfg$nrep),
            file.path(results_dir,"manifests","runtime_summary.csv"))
  list(results_dir=results_dir, cell_performance=cell_perf)
}
