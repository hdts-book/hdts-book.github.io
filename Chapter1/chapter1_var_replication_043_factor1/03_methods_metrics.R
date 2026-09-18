# sparseVAR 0.4.3 method family and performance measures.
# Production profile: standardize=TRUE, penalize_diag=FALSE, whiten='none'.

forecast_errors <- function(Phi, Y_sample, Y_test) {
  if (any(!is.finite(Phi))) return(list(msfe = NA_real_, rmfe = NA_real_))
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
  f1den <- 2 * TP + FP + FN
  F1 <- if (f1den > 0) 2 * TP / f1den else NA_real_
  list(TP=TP, FP=FP, FN=FN, TN=TN, precision=precision, recall=recall, F1=F1)
}

metrics_one <- function(Phi, Atrue, Y_sample, Y_test, cfg) {
  d <- nrow(Atrue)
  if (is.null(Phi) || is.null(dim(Phi)) || !all(dim(Phi) == c(d, d)) || any(!is.finite(Phi))) {
    stop("Invalid/non-finite coefficient matrix returned by method.")
  }

  S0 <- abs(Atrue) > cfg$support_tol
  Sh <- abs(Phi) > cfg$support_tol
  cc <- support_counts(S0, Sh)
  specificity <- if (cc$TN + cc$FP > 0) cc$TN / (cc$TN + cc$FP) else NA_real_

  off <- row(Atrue) != col(Atrue)
  S0o <- S0[off]; Sho <- Sh[off]
  oo <- support_counts(S0o, Sho)
  # F1/recall are not meaningful when the true off-diagonal support is empty.
  # Keep off-diagonal Hamming/FP diagnostics, but mark these two as NA so an
  # empty-support DGP cell cannot distort off-diagonal F1 ranks.
  if (sum(S0o) == 0L) {
    oo$recall <- NA_real_
    oo$F1 <- NA_real_
  }
  off_n <- d * (d - 1L)

  fcast <- forecast_errors(Phi, Y_sample, Y_test)
  fsq <- sum((Phi - Atrue)^2)

  data.frame(
    frob = sqrt(fsq),
    frob_sq = fsq,
    msfe = fcast$msfe,
    rmfe = fcast$rmfe,
    hamming = (cc$FP + cc$FN) / (d * d),
    precision = cc$precision,
    recall = cc$recall,
    F1 = cc$F1,
    specificity = specificity,
    TPR = cc$recall,
    FPR = if (cc$FP + cc$TN > 0) cc$FP / (cc$FP + cc$TN) else NA_real_,
    TP = cc$TP, FP = cc$FP, FN = cc$FN, TN = cc$TN,
    support_true = sum(S0),
    support_hat = sum(Sh),
    off_hamming = if (off_n > 0) (oo$FP + oo$FN) / off_n else NA_real_,
    off_precision = oo$precision,
    off_recall = oo$recall,
    off_F1 = oo$F1,
    off_TP = oo$TP, off_FP = oo$FP, off_FN = oo$FN, off_TN = oo$TN,
    support_true_offdiag = sum(S0o),
    support_hat_offdiag = sum(Sho),
    rho_true = spectral_radius(Atrue),
    rho_hat = spectral_radius(Phi),
    max_abs = max(abs(Phi)),
    stringsAsFactors = FALSE
  )
}

fit_method_family_043 <- function(Y_sample, tuning, cfg, method_order) {
  tuning <- match.arg(tuning, cfg$tunings)
  t0 <- proc.time()[3L]
  fit <- sparseVAR::sVAR_variants(
    Yt = Y_sample,
    p = cfg$p,
    tuning = tuning,
    fold = cfg$fold,
    theory_c = cfg$theory_c,
    penalize_diag = cfg$penalize_diag,
    standardize = cfg$standardize,
    whiten = cfg$whiten,
    threshold_factor = cfg$threshold_factor,
    relax_gamma = cfg$relax_gamma,
    scaled_lambda0 = cfg$scaled_lambda0,
    nlambda = cfg$nlambda,
    max_iter = cfg$max_iter,
    tol = cfg$tol,
    adaptive_zero = cfg$adaptive_zero,
    adaptive_cv = cfg$adaptive_cv
  )
  runtime <- proc.time()[3L] - t0
  assert_method_set(fit$methods, method_order)
  fit$methods <- fit$methods[method_order]

  if (!identical(fit$base_fit$whiten, cfg$whiten)) stop("Package returned an unexpected whitening profile.")
  if (!identical(fit$base_fit$penalize_diag, cfg$penalize_diag)) stop("Package returned an unexpected diagonal-penalty profile.")
  if (!is.null(fit$base_fit$SigmaInv)) stop("whiten='none' must return SigmaInv=NULL.")
  list(fit = fit, runtime = runtime)
}

method_replication_043 <- function(
    dgp, idx_d, idx_T, rr, Atrue, data_seed_value,
    cfg, method_order) {

  d <- cfg$d_grid[idx_d]
  TT <- cfg$T_grid[idx_T]
  dat <- simulate_var1(Atrue, TT, data_seed_value, cfg)
  truth_seed_value <- truth_seed(dgp, idx_d, cfg$truth_seed_base)

  out <- vector("list", length(cfg$tunings))
  for (tt in seq_along(cfg$tunings)) {
    tuning <- cfg$tunings[tt]
    obj <- fit_method_family_043(dat$sample, tuning, cfg, method_order)
    fit <- obj$fit
    base <- fit$base_fit

    lasso_support <- sum(abs(base$Phi_hat_lasso) > cfg$support_tol)
    adaptive_support <- sum(abs(base$Phi_hat_Alasso) > cfg$support_tol)
    off <- row(Atrue) != col(Atrue)
    lasso_support_off <- sum(abs(base$Phi_hat_lasso[off]) > cfg$support_tol)
    adaptive_support_off <- sum(abs(base$Phi_hat_Alasso[off]) > cfg$support_tol)
    scaled_lambda0 <- sqrt(2 * log(d * cfg$p) / (TT - cfg$p))

    rows <- do.call(rbind, lapply(method_order, function(mm) {
      m <- metrics_one(fit$methods[[mm]], Atrue, dat$sample, dat$test, cfg)
      m$rep <- rr
      m$dgp <- dgp
      m$idx_d <- idx_d
      m$d <- d
      m$idx_T <- idx_T
      m$T_sample <- TT
      m$tuning <- tuning
      m$method <- mm
      m$data_seed <- data_seed_value
      m$truth_seed <- truth_seed_value
      m$lambda_lasso <- base$lambda_lasso
      m$lambda_adaptive <- if (startsWith(mm, "A")) base$lambda_adalasso else NA_real_
      m$scaled_lambda0 <- if (mm == "ScL") scaled_lambda0 else NA_real_
      m$base_lasso_support <- lasso_support
      m$base_adaptive_support <- adaptive_support
      m$base_lasso_support_offdiag <- lasso_support_off
      m$base_adaptive_support_offdiag <- adaptive_support_off
      m$threshold_factor <- cfg$threshold_factor
      m$runtime_fit_seconds <- obj$runtime
      m
    }))

    rows <- rows[, c(
      "rep", "dgp", "idx_d", "d", "idx_T", "T_sample", "tuning", "method",
      "data_seed", "truth_seed",
      "frob", "frob_sq", "msfe", "rmfe", "hamming", "precision", "recall", "F1",
      "specificity", "TPR", "FPR", "TP", "FP", "FN", "TN",
      "support_true", "support_hat",
      "off_hamming", "off_precision", "off_recall", "off_F1",
      "off_TP", "off_FP", "off_FN", "off_TN",
      "support_true_offdiag", "support_hat_offdiag",
      "rho_true", "rho_hat", "max_abs",
      "lambda_lasso", "lambda_adaptive", "scaled_lambda0",
      "base_lasso_support", "base_adaptive_support",
      "base_lasso_support_offdiag", "base_adaptive_support_offdiag",
      "threshold_factor", "runtime_fit_seconds"
    )]

    core <- c("frob", "rmfe", "hamming")
    if (any(!is.finite(as.matrix(rows[, core, drop = FALSE])))) {
      stop("Non-finite core metric: ", dgp, " d=", d, " T=", TT,
           " rep=", rr, " tuning=", tuning)
    }
    if (any(is.na(rows$F1))) {
      stop("Undefined full-matrix F1 encountered: ", dgp, " d=", d, " T=", TT,
           " rep=", rr, " tuning=", tuning)
    }
    out[[tt]] <- rows
  }

  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}
