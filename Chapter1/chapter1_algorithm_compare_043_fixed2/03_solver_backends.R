# Common objective and solver wrappers.
#
# In standardized coordinates all solvers target
#   (1/N) ||Y - X B||_F^2 + lambda * sum_ij W_ij |B_ij|,
# with W_ii=0 for VAR(1) own-lag diagonal coefficients.
# No residual-covariance whitening is used.

prepare_common_regression <- function(Yt, cfg) {
  Yc <- as.matrix(Yt)
  Yc <- Yc - rowMeans(Yc)
  des <- getFromNamespace("var_design", "sparseVAR")(Yc, p = cfg$p)
  X <- as.matrix(des$X_design)
  Y <- as.matrix(des$Y_resp)
  N <- nrow(X); d <- ncol(X); k <- ncol(Y)
  X_sd <- rep(1, d); Y_sd <- rep(1, k)
  Xs <- X; Ys <- Y
  if (isTRUE(cfg$standardize)) {
    X_sd <- apply(X, 2, stats::sd)
    Y_sd <- apply(Y, 2, stats::sd)
    X_sd[!is.finite(X_sd) | X_sd <= 0] <- 1
    Y_sd[!is.finite(Y_sd) | Y_sd <= 0] <- 1
    Xs <- sweep(X, 2, X_sd, "/")
    Ys <- sweep(Y, 2, Y_sd, "/")
  }
  list(Y_centered = Yc, X = X, Y = Y, Xs = as.matrix(Xs), Ys = as.matrix(Ys),
       X_sd = X_sd, Y_sd = Y_sd, N = N, d = d, k = k)
}

make_penalty_weights <- function(d, k, penalize_diag = FALSE) {
  W <- matrix(1, d, k)
  if (!isTRUE(penalize_diag)) {
    for (j in seq_len(min(d, k))) W[j, j] <- 0
  }
  W
}

backtransform_B <- function(B, pre) {
  diag(pre$Y_sd, pre$k) %*% t(B) %*% diag(1 / pre$X_sd, pre$d)
}

book_theory_lambda <- function(pre, cfg) {
  q <- cfg$p * pre$k^2
  cfg$theory_c * sqrt(log(max(q, 2)) / pre$N)
}

fit_fista_B <- function(X, Y, W, lambda, cfg, pilot = FALSE) {
  fn <- getFromNamespace("fista_lasso_multi_cpp", "sparseVAR")
  mi <- if (pilot) cfg$pilot_max_iter else cfg$fista_max_iter
  tt <- if (pilot) cfg$pilot_tol else cfg$fista_tol
  t0 <- proc.time()[3L]
  B <- fn(as.matrix(X), as.matrix(Y), as.matrix(W), as.numeric(lambda),
          matrix(0, 0, 0), as.integer(mi), as.numeric(tt))
  list(B = as.matrix(B), runtime = proc.time()[3L] - t0,
       converged = NA, n_iter = NA_integer_)
}

fit_admm_B <- function(X, Y, W, lambda, cfg) {
  fn <- getFromNamespace("admm_varp_multi_cpp", "sparseVAR")
  t0 <- proc.time()[3L]
  out <- fn(as.matrix(X), as.matrix(Y), as.matrix(W), as.numeric(lambda),
            matrix(0, 0, 0), as.integer(cfg$admm_max_iter),
            as.numeric(cfg$admm_rho), as.numeric(cfg$admm_tol))
  list(B = as.matrix(out$B), runtime = proc.time()[3L] - t0,
       converged = isTRUE(out$converged), n_iter = as.integer(out$n_iter))
}

# glmnet solves (1/(2N))RSS + lambda_g * sum_j pf'_j |beta_j| and internally
# rescales penalty.factor to sum to the number of predictors.  We fit equations
# separately.  The mapping below makes its effective penalty exactly
# (lambda/2) * W[,j], which is the common book objective after multiplying by 1/2.
fit_glmnet_B <- function(X, Y, W, lambda, cfg) {
  X <- as.matrix(X); Y <- as.matrix(Y); W <- as.matrix(W)
  d <- ncol(X); k <- ncol(Y)
  B <- matrix(0, d, k)
  t0 <- proc.time()[3L]
  for (j in seq_len(k)) {
    pf <- as.numeric(W[, j])
    sw <- sum(pf)
    if (!is.finite(sw) || sw <= 0) {
      B[, j] <- as.numeric(solve(crossprod(X) + 1e-12 * diag(d), crossprod(X, Y[, j])))
      next
    }
    lambda_g <- (lambda / 2) * (sw / d)
    fit <- glmnet::glmnet(
      x = X, y = as.numeric(Y[, j]), family = "gaussian",
      alpha = 1, lambda = lambda_g,
      intercept = FALSE, standardize = FALSE,
      penalty.factor = pf,
      thresh = cfg$glmnet_thresh, maxit = cfg$glmnet_maxit
    )
    B[, j] <- as.numeric(stats::coef(fit, s = lambda_g))[-1L]
  }
  list(B = B, runtime = proc.time()[3L] - t0,
       converged = NA, n_iter = NA_integer_)
}


canonicalize_solver_B <- function(B, cfg) {
  B <- as.matrix(B)
  idx <- abs(B) <= cfg$solver_zero_tol
  cleaned <- abs(B[idx])
  max_cleaned <- if (length(cleaned)) max(cleaned) else 0
  n_cleaned <- sum(idx & (B != 0))
  B[idx] <- 0
  list(B = B, n_cleaned = n_cleaned, max_cleaned_abs = max_cleaned)
}

fit_solver_B <- function(solver, X, Y, W, lambda, cfg) {
  out <- switch(solver,
                glmnet = fit_glmnet_B(X, Y, W, lambda, cfg),
                admm = fit_admm_B(X, Y, W, lambda, cfg),
                fista = fit_fista_B(X, Y, W, lambda, cfg, pilot = FALSE),
                stop("Unknown solver: ", solver))
  cc <- canonicalize_solver_B(out$B, cfg)
  out$B_raw <- out$B
  out$B <- cc$B
  out$n_cleaned <- cc$n_cleaned
  out$max_cleaned_abs <- cc$max_cleaned_abs
  out
}

adaptive_weights_common <- function(Bpilot, W0, cfg) {
  W <- 1 / pmax(abs(Bpilot), cfg$adaptive_eps)
  W <- pmin(W, cfg$adaptive_cap)
  W[W0 == 0] <- 0
  W
}

reference_tuning <- function(Y_sample, tuning, pre, W0, cfg) {
  if (tuning == "theory") {
    lam <- book_theory_lambda(pre, cfg)
    return(list(lambda_lasso = lam, lambda_adaptive = lam,
                tuning_runtime = 0, source = "theory"))
  }
  t0 <- proc.time()[3L]
  fit <- sparseVAR::sVAR_adalasso_fista(
    Yt = Y_sample, p = cfg$p, fold = cfg$fold,
    penalize_diag = cfg$penalize_diag,
    standardize = cfg$standardize,
    whiten = cfg$whiten,
    lambda_rule = "cv", theory_c = cfg$theory_c,
    nlambda = cfg$nlambda, lambda_min_ratio = cfg$lambda_min_ratio,
    max_iter = cfg$reference_max_iter, tol = cfg$reference_tol,
    adaptive_zero = cfg$adaptive_zero, adaptive_cv = cfg$adaptive_cv
  )
  if (!is.null(fit$SigmaInv)) stop("Reference CV unexpectedly used whitening.")
  list(lambda_lasso = fit$lambda_lasso,
       lambda_adaptive = fit$lambda_adalasso,
       tuning_runtime = proc.time()[3L] - t0,
       source = "fista_nested_cv_reference")
}

objective_book <- function(B, X, Y, W, lambda) {
  R <- Y - X %*% B
  sum(R^2) / nrow(X) + lambda * sum(W * abs(B))
}

kkt_violation <- function(B, X, Y, W, lambda, zero_tol = 1e-6) {
  G <- (2 / nrow(X)) * crossprod(X, X %*% B - Y)
  out <- matrix(0, nrow(B), ncol(B))
  active <- abs(B) > zero_tol
  pen <- W > 0
  idx <- active & pen
  out[idx] <- abs(G[idx] + lambda * W[idx] * sign(B[idx]))
  idx0 <- (!active) & pen
  out[idx0] <- pmax(abs(G[idx0]) - lambda * W[idx0], 0)
  unpen <- !pen
  out[unpen] <- abs(G[unpen])
  max(out)
}
