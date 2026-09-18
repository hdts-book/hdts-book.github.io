# Statistical metrics and one-replication solver comparison.

forecast_errors <- function(Phi, Y_sample, Y_test) {
  Yall <- cbind(Y_sample, Y_test)
  T0 <- ncol(Y_sample)
  e2 <- numeric(ncol(Y_test))
  for (h in seq_len(ncol(Y_test))) {
    e <- Y_test[, h] - Phi %*% Yall[, T0 + h - 1L]
    e2[h] <- sum(e^2)
  }
  list(msfe = mean(e2), rmfe = sqrt(mean(e2)))
}

support_counts <- function(S0, Sh) {
  TP <- sum(Sh & S0); FP <- sum(Sh & !S0)
  FN <- sum(!Sh & S0); TN <- sum(!Sh & !S0)
  precision <- if (TP + FP > 0) TP / (TP + FP) else 0
  recall <- if (TP + FN > 0) TP / (TP + FN) else 0
  den <- 2 * TP + FP + FN
  F1 <- if (den > 0) 2 * TP / den else NA_real_
  list(TP=TP, FP=FP, FN=FN, TN=TN, precision=precision, recall=recall, F1=F1)
}

metrics_one <- function(Phi, Atrue, Y_sample, Y_test, cfg) {
  d <- nrow(Atrue)
  S0 <- abs(Atrue) > cfg$support_tol
  Sh <- abs(Phi) > cfg$support_tol
  cc <- support_counts(S0, Sh)
  off <- row(Atrue) != col(Atrue)
  oo <- support_counts(S0[off], Sh[off])
  if (sum(S0[off]) == 0L) { oo$recall <- NA_real_; oo$F1 <- NA_real_ }
  fsq <- sum((Phi - Atrue)^2)
  fc <- forecast_errors(Phi, Y_sample, Y_test)
  data.frame(
    frob = sqrt(fsq), frob_sq = fsq, msfe = fc$msfe, rmfe = fc$rmfe,
    hamming = (cc$FP + cc$FN) / (d*d),
    precision = cc$precision, recall = cc$recall, F1 = cc$F1,
    TP = cc$TP, FP = cc$FP, FN = cc$FN, TN = cc$TN,
    support_true = sum(S0), support_hat = sum(Sh),
    off_hamming = (oo$FP + oo$FN) / (d*(d-1L)),
    off_precision = oo$precision, off_recall = oo$recall, off_F1 = oo$F1,
    off_TP = oo$TP, off_FP = oo$FP, off_FN = oo$FN, off_TN = oo$TN,
    support_true_offdiag = sum(S0[off]), support_hat_offdiag = sum(Sh[off]),
    rho_true = spectral_radius(Atrue), rho_hat = spectral_radius(Phi),
    stringsAsFactors = FALSE
  )
}

one_algorithm_replication <- function(dgp, idx_d, idx_T, rr, Atrue,
                                      data_seed_value, cfg, solver_order, family_order) {
  d <- cfg$d_grid[idx_d]; TT <- cfg$T_grid[idx_T]
  dat <- simulate_var1(Atrue, TT, data_seed_value, cfg)
  pre <- prepare_common_regression(dat$sample, cfg)
  W0 <- make_penalty_weights(pre$d, pre$k, cfg$penalize_diag)
  truth_seed_value <- truth_seed(dgp, idx_d, cfg$truth_seed_base)

  all_tuning <- vector("list", length(cfg$tunings))
  for (it in seq_along(cfg$tunings)) {
    tuning <- cfg$tunings[it]
    ref <- reference_tuning(dat$sample, tuning, pre, W0, cfg)

    # Common high-accuracy pilot. This is intentionally outside the solver
    # comparison so adaptive weights are identical for glmnet/ADMM/FISTA.
    p0 <- fit_fista_B(pre$Xs, pre$Ys, W0, ref$lambda_lasso, cfg, pilot = TRUE)
    W1 <- adaptive_weights_common(p0$B, W0, cfg)

    family_rows <- vector("list", length(family_order))
    for (ifam in seq_along(family_order)) {
      fam <- family_order[ifam]
      W <- if (fam == "lasso") W0 else W1
      lam <- if (fam == "lasso") ref$lambda_lasso else ref$lambda_adaptive

      fits <- lapply(solver_order, function(s) fit_solver_B(s, pre$Xs, pre$Ys, W, lam, cfg))
      names(fits) <- solver_order
      Bf <- fits[["fista"]]$B
      obj <- vapply(fits, function(z) objective_book(z$B, pre$Xs, pre$Ys, W, lam), numeric(1))
      best_obj <- min(obj)

      rows <- do.call(rbind, lapply(solver_order, function(s) {
        z <- fits[[s]]
        Phi <- backtransform_B(z$B, pre)
        m <- metrics_one(Phi, Atrue, dat$sample, dat$test, cfg)
        m$rep <- rr; m$dgp <- dgp; m$idx_d <- idx_d; m$d <- d
        m$idx_T <- idx_T; m$T_sample <- TT
        m$tuning <- tuning; m$family <- fam; m$solver <- s
        m$method <- paste(fam, s, sep = "_")
        m$data_seed <- data_seed_value; m$truth_seed <- truth_seed_value
        m$lambda_lasso <- ref$lambda_lasso
        m$lambda_adaptive <- ref$lambda_adaptive
        m$lambda_used <- lam
        m$reference_tuning_source <- ref$source
        m$reference_tuning_seconds <- ref$tuning_runtime
        m$reference_pilot_seconds <- p0$runtime
        m$runtime_seconds <- z$runtime
        m$objective <- obj[[s]]
        m$objective_gap <- obj[[s]] - best_obj
        m$kkt_max <- kkt_violation(z$B, pre$Xs, pre$Ys, W, lam, zero_tol = cfg$solver_zero_tol)
        m$n_solver_zeroed <- z$n_cleaned
        m$max_solver_zeroed_abs <- z$max_cleaned_abs
        m$B_diff_to_fista <- sqrt(sum((z$B - Bf)^2))
        m$Phi_diff_to_fista <- sqrt(sum((Phi - backtransform_B(Bf, pre))^2))
        m$admm_converged <- if (s == "admm") isTRUE(z$converged) else NA
        m$solver_n_iter <- if (s == "admm") z$n_iter else NA_integer_
        m
      }))
      family_rows[[ifam]] <- rows
    }
    all_tuning[[it]] <- do.call(rbind, family_rows)
  }

  ans <- do.call(rbind, all_tuning)
  core <- c("frob", "rmfe", "hamming", "F1", "objective", "kkt_max")
  if (any(!is.finite(as.matrix(ans[, core, drop = FALSE])))) {
    stop("Non-finite solver-comparison metric: ", dgp, " d=", d, " T=", TT, " rep=", rr)
  }
  rownames(ans) <- NULL
  ans
}
