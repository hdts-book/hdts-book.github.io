# Utilities for the independent-design VR replication.
# sparseVAR 0.4.3 supplies the FISTA backend; the VR experiment itself is an
# independent multivariate regression benchmark rather than a time-series VAR fit.

`%||%` <- function(x, y) if (is.null(x)) y else x

set_rng_seed <- function(seed, cfg = CFG) {
  RNGversion(cfg$rng_version)
  RNGkind(
    kind = cfg$rng_kind,
    normal.kind = cfg$rng_normal_kind,
    sample.kind = cfg$rng_sample_kind
  )
  set.seed(as.integer(seed))
  invisible(seed)
}

largest_singular_value <- function(A) {
  svd(A, nu = 0L, nv = 0L)$d[1L]
}

spectral_radius <- function(A) {
  if (any(!is.finite(A))) return(NA_real_)
  max(Mod(eigen(A, only.values = TRUE)$values))
}

truth_seed <- function(dgp, idx_d, base) {
  off <- switch(
    dgp,
    dgp1 = 0L,
    dgp2 = 1000L,
    dgp3 = 2000L,
    stop("Unknown DGP: ", dgp)
  )
  as.integer(base + off + idx_d)
}

data_seed <- function(master_seed, dgp, idx_d, idx_T, rr) {
  off <- switch(
    dgp,
    dgp1 = 0L,
    dgp2 = 1000000L,
    dgp3 = 2000000L,
    stop("Unknown DGP: ", dgp)
  )
  as.integer(master_seed + off + idx_d * 100000L + idx_T * 1000L + rr)
}

make_block_sizes <- function(d, K = 2L) {
  z <- rep(floor(d / K), K)
  z[K] <- d - sum(z[-K])
  z
}

make_A_true <- function(dgp, d, idx_d, cfg) {
  set_rng_seed(truth_seed(dgp, idx_d, cfg$truth_seed_base), cfg)

  if (dgp == "dgp1") {
    M <- matrix(0, d, d)
    oi <- which(row(M) != col(M))
    M[oi] <- stats::rbinom(length(oi), 1L, cfg$dgp1_prop)
    diag(M) <- 1
  } else if (dgp == "dgp2") {
    M <- matrix(0, d, d)
    diag(M) <- stats::runif(d, cfg$dgp2_a1, cfg$dgp2_a2)
    if (d >= 2L) {
      for (i in seq_len(d - 1L)) {
        M[i, i + 1L] <- stats::runif(1L, cfg$dgp2_b1, cfg$dgp2_b2)
        M[i + 1L, i] <- stats::runif(1L, cfg$dgp2_b1, cfg$dgp2_b2)
      }
    }
  } else if (dgp == "dgp3") {
    bs <- make_block_sizes(d, cfg$dgp3_K)
    P <- matrix(cfg$dgp3_p_out, cfg$dgp3_K, cfg$dgp3_K)
    diag(P) <- cfg$dgp3_p_in
    g <- igraph::sample_sbm(
      n = d,
      pref.matrix = P,
      block.sizes = bs,
      directed = TRUE,
      loops = FALSE
    )
    # Historical Chapter-1 orientation.
    M <- t(as.matrix(igraph::as_adj(g, sparse = FALSE)))
    diag(M) <- 1
  } else {
    stop("Unknown DGP: ", dgp)
  }

  sm <- largest_singular_value(M)
  if (!is.finite(sm) || sm <= 0) stop("Degenerate DGP matrix.")
  cfg$alpha * M / sm
}

stationary_cov_var1 <- function(A, Sigma_u = NULL) {
  d <- nrow(A)
  if (is.null(Sigma_u)) Sigma_u <- diag(d)
  K <- diag(d * d) - kronecker(A, A)
  g <- solve(K, as.vector(Sigma_u))
  G <- matrix(g, d, d)
  (G + t(G)) / 2
}

signal_scale_for_snr <- function(A_base, Sigma_X, Sigma_e, snr_target) {
  signal <- sum(diag(A_base %*% Sigma_X %*% t(A_base)))
  noise <- sum(diag(Sigma_e))
  if (!is.finite(signal) || signal <= 0) stop("Invalid signal variance.")
  sqrt(snr_target * noise / signal)
}

simulate_vr <- function(A_base, T_sample, T_test, snr_target, seed, cfg) {
  set_rng_seed(seed, cfg)
  d <- nrow(A_base)
  Sigma_e <- diag(d)

  # Historical VR design: X has the stationary covariance induced by the VAR
  # benchmark, but rows are independent over observations.
  Sigma_X <- stationary_cov_var1(A_base, Sigma_e)
  c_sig <- signal_scale_for_snr(A_base, Sigma_X, Sigma_e, snr_target)
  A_true <- c_sig * A_base

  n <- T_sample + T_test
  X <- MASS::mvrnorm(n, mu = rep(0, d), Sigma = Sigma_X)
  E <- MASS::mvrnorm(n, mu = rep(0, d), Sigma = Sigma_e)
  Y <- X %*% t(A_true) + E

  list(
    X_sample = X[seq_len(T_sample), , drop = FALSE],
    Y_sample = Y[seq_len(T_sample), , drop = FALSE],
    X_test = X[T_sample + seq_len(T_test), , drop = FALSE],
    Y_test = Y[T_sample + seq_len(T_test), , drop = FALSE],
    A_true = A_true,
    A_base = A_base,
    signal_scale = c_sig,
    Sigma_X = Sigma_X,
    Sigma_e = Sigma_e
  )
}

block_foldid <- function(n, k) {
  floor((seq_len(n) - 1L) * k / n) + 1L
}

safe_invsqrt <- function(S) {
  ee <- eigen((S + t(S)) / 2, symmetric = TRUE)
  vv <- pmax(ee$values, 1e-8)
  ee$vectors %*% diag(1 / sqrt(vv), length(vv)) %*% t(ee$vectors)
}

prepare_vr <- function(X, Y, cfg) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  N <- nrow(X)
  d <- ncol(X)
  k <- ncol(Y)
  if (nrow(Y) != N) stop("X and Y must have the same number of rows.")

  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- scale(Y, center = TRUE, scale = FALSE)

  X_sd <- rep(1, d)
  Y_sd <- rep(1, k)
  Xs <- Xc
  Ys <- Yc

  if (isTRUE(cfg$standardize)) {
    X_sd <- apply(Xc, 2, stats::sd)
    Y_sd <- apply(Yc, 2, stats::sd)
    X_sd[!is.finite(X_sd) | X_sd <= 0] <- 1
    Y_sd[!is.finite(Y_sd) | Y_sd <= 0] <- 1
    Xs <- sweep(Xc, 2, X_sd, "/")
    Ys <- sweep(Yc, 2, Y_sd, "/")
  }

  Sigma_hat <- NULL
  SigmaInv <- NULL

  if (identical(cfg$whiten, "none")) {
    if (isTRUE(cfg$updateSigma)) {
      stop("whiten='none' must be paired with updateSigma=FALSE.")
    }
  } else {
    if (!isTRUE(cfg$updateSigma)) {
      stop("whiten='diag'/'full' requires updateSigma=TRUE.")
    }
    B0 <- solve(crossprod(Xs) + 1e-8 * diag(d), crossprod(Xs, Ys))
    R0 <- Ys - Xs %*% B0
    Sigma_hat <- crossprod(R0) / N
    if (identical(cfg$whiten, "diag")) {
      SigmaInv <- diag(1 / pmax(diag(Sigma_hat), 1e-8), k)
    } else if (identical(cfg$whiten, "full")) {
      H <- safe_invsqrt(Sigma_hat)
      SigmaInv <- t(H) %*% H
    } else {
      stop("Unknown whitening mode: ", cfg$whiten)
    }
  }

  list(
    X = as.matrix(Xc),
    Y = as.matrix(Yc),
    Xs = as.matrix(Xs),
    Ys = as.matrix(Ys),
    X_sd = as.numeric(X_sd),
    Y_sd = as.numeric(Y_sd),
    Sigma_hat = Sigma_hat,
    SigmaInv = SigmaInv,
    N = N,
    d = d,
    k = k
  )
}

fista_backend <- function(X, Y, W, lambda, SigmaInv, cfg) {
  fn <- getFromNamespace("fista_lasso_multi_cpp", "sparseVAR")
  Sarg <- if (is.null(SigmaInv)) matrix(0, 0, 0) else as.matrix(SigmaInv)
  out <- fn(
    as.matrix(X), as.matrix(Y), as.matrix(W), as.numeric(lambda),
    Sarg, as.integer(cfg$max_iter), as.numeric(cfg$tol)
  )
  if (is.list(out) && !is.null(out$B)) out <- out$B
  out <- as.matrix(out)
  if (any(!is.finite(out))) stop("Non-finite FISTA coefficient matrix.")
  out
}

weighted_sse <- function(R, SigmaInv) {
  if (is.null(SigmaInv)) sum(R^2) else sum((R %*% SigmaInv) * R)
}

make_penalty_weights <- function(d, k, penalize_diag) {
  W <- matrix(1, d, k)
  if (!isTRUE(penalize_diag)) {
    for (j in seq_len(min(d, k))) W[j, j] <- 0
  }
  W
}

kkt_lambda_max <- function(X, Y, W, SigmaInv) {
  G <- 2 * crossprod(X, Y) / nrow(X)
  if (!is.null(SigmaInv)) G <- G %*% SigmaInv
  den <- W
  den[den <= 0] <- Inf
  z <- abs(G) / den
  z[!is.finite(z)] <- NA_real_
  ans <- suppressWarnings(max(z, na.rm = TRUE))
  if (!is.finite(ans) || ans <= 0) ans <- 1
  ans
}

lambda_grid <- function(X, Y, W, SigmaInv, cfg) {
  lm <- kkt_lambda_max(X, Y, W, SigmaInv)
  exp(seq(
    log(lm),
    log(lm * cfg$lambda_min_ratio),
    length.out = cfg$nlambda
  ))
}

cv_lambda <- function(X, Y, W, grid, SigmaInv, cfg) {
  fid <- block_foldid(nrow(X), cfg$fold)
  loss <- numeric(length(grid))

  for (f in seq_len(cfg$fold)) {
    tr <- which(fid != f)
    va <- which(fid == f)
    for (ii in seq_along(grid)) {
      B <- fista_backend(
        X[tr, , drop = FALSE],
        Y[tr, , drop = FALSE],
        W, grid[ii], SigmaInv, cfg
      )
      R <- Y[va, , drop = FALSE] - X[va, , drop = FALSE] %*% B
      loss[ii] <- loss[ii] + weighted_sse(R, SigmaInv)
    }
  }
  grid[which.min(loss)]
}

adaptive_weights <- function(B, W0, cfg) {
  W <- 1 / pmax(abs(B), cfg$adaptive_eps)
  W <- pmin(W, cfg$adaptive_cap)
  # Structural zero penalty entries (diagonal when penalize_diag=FALSE) must
  # remain unpenalized in the adaptive stage as well.
  W[W0 == 0] <- 0
  W
}

nested_adaptive_lambda <- function(X, Y, W0, grid2, SigmaInv, cfg) {
  outer <- block_foldid(nrow(X), cfg$fold)
  loss <- numeric(length(grid2))

  for (f in seq_len(cfg$fold)) {
    tr <- which(outer != f)
    va <- which(outer == f)

    grid1 <- lambda_grid(
      X[tr, , drop = FALSE],
      Y[tr, , drop = FALSE],
      W0, SigmaInv, cfg
    )
    lambda1 <- cv_lambda(
      X[tr, , drop = FALSE],
      Y[tr, , drop = FALSE],
      W0, grid1, SigmaInv, cfg
    )
    B1 <- fista_backend(
      X[tr, , drop = FALSE],
      Y[tr, , drop = FALSE],
      W0, lambda1, SigmaInv, cfg
    )
    W1 <- adaptive_weights(B1, W0, cfg)

    for (ii in seq_along(grid2)) {
      B2 <- fista_backend(
        X[tr, , drop = FALSE],
        Y[tr, , drop = FALSE],
        W1, grid2[ii], SigmaInv, cfg
      )
      R <- Y[va, , drop = FALSE] - X[va, , drop = FALSE] %*% B2
      loss[ii] <- loss[ii] + weighted_sse(R, SigmaInv)
    }
  }
  grid2[which.min(loss)]
}

backtransform <- function(B, pre) {
  diag(pre$Y_sd, pre$k) %*% t(B) %*% diag(1 / pre$X_sd, pre$d)
}

fit_vr_fista <- function(X, Y, tuning, cfg) {
  tuning <- match.arg(tuning, cfg$tunings)
  pre <- prepare_vr(X, Y, cfg)
  q <- pre$d * pre$k
  W0 <- make_penalty_weights(pre$d, pre$k, cfg$penalize_diag)

  grid1 <- lambda_grid(pre$Xs, pre$Ys, W0, pre$SigmaInv, cfg)
  if (tuning == "theory") {
    lambda1 <- cfg$theory_c * sqrt(log(q) / pre$N)
  } else {
    lambda1 <- cv_lambda(pre$Xs, pre$Ys, W0, grid1, pre$SigmaInv, cfg)
  }

  B1 <- fista_backend(pre$Xs, pre$Ys, W0, lambda1, pre$SigmaInv, cfg)
  W1 <- adaptive_weights(B1, W0, cfg)
  grid2 <- lambda_grid(pre$Xs, pre$Ys, W1, pre$SigmaInv, cfg)

  if (tuning == "theory") {
    if (!isTRUE(cfg$theory_same_lambda)) {
      stop("This production code expects theory_same_lambda=TRUE.")
    }
    lambda2 <- lambda1
  } else {
    lambda2 <- nested_adaptive_lambda(
      pre$Xs, pre$Ys, W0, grid2, pre$SigmaInv, cfg
    )
  }

  B2 <- fista_backend(pre$Xs, pre$Ys, W1, lambda2, pre$SigmaInv, cfg)

  list(
    Phi_hat_lasso = backtransform(B1, pre),
    Phi_hat_Alasso = backtransform(B2, pre),
    lambda_lasso = lambda1,
    lambda_adalasso = lambda2,
    pre = pre,
    W0 = W0,
    W1 = W1,
    B_lasso = B1,
    B_adaptive = B2
  )
}

vr_refit <- function(X, Y, Phi) {
  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  d <- ncol(Xc)
  k <- ncol(Yc)
  out <- matrix(0, k, d)

  for (j in seq_len(k)) {
    idx <- which(abs(Phi[j, ]) > 0)
    if (!length(idx)) next
    XtX <- crossprod(Xc[, idx, drop = FALSE])
    rhs <- crossprod(Xc[, idx, drop = FALSE], Yc[, j])
    b <- solve(XtX + 1e-8 * diag(length(idx)), rhs)
    out[j, idx] <- as.numeric(b)
  }
  out
}

vr_debias <- function(X, Y, Phi) {
  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  N <- nrow(Xc)
  B <- t(Phi)
  G <- crossprod(Xc) / N
  M <- solve(G + 1e-10 * diag(ncol(G)))
  score <- crossprod(Xc, Yc - Xc %*% B) / N
  t(B + M %*% score)
}

vr_threshold <- function(Phi, eta, cfg) {
  out <- as.matrix(Phi)
  d <- nrow(out)
  k <- ncol(out)
  mask <- abs(out) < eta
  if (!isTRUE(cfg$penalize_diag)) {
    for (j in seq_len(min(d, k))) mask[j, j] <- FALSE
  }
  out[mask] <- 0
  out
}

vr_scaled_lasso <- function(X, Y, cfg) {
  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  N <- nrow(Xc)
  d <- ncol(Xc)
  k <- ncol(Yc)

  xsd <- apply(Xc, 2, stats::sd)
  xsd[!is.finite(xsd) | xsd <= 0] <- 1
  Xs <- if (isTRUE(cfg$standardize)) sweep(Xc, 2, xsd, "/") else Xc
  if (!isTRUE(cfg$standardize)) xsd[] <- 1

  lambda0 <- cfg$scaled_lambda0
  if (is.null(lambda0)) lambda0 <- cfg$theory_c * sqrt(2 * log(d) / N)

  B <- matrix(0, d, k)
  sigma <- sqrt(colMeans(Yc^2))
  sigma[!is.finite(sigma) | sigma <= 1e-10] <- 1
  nit <- integer(k)

  for (j in seq_len(k)) {
    # Under standardization, coefficients in the scaled design receive weights
    # 1/xsd so the penalty corresponds to the original coefficient scale.
    w <- matrix(1 / xsd, d, 1)
    if (!isTRUE(cfg$penalize_diag) && j <= d) w[j, 1] <- 0
    sj <- sigma[j]
    beta <- rep(0, d)

    for (it in seq_len(cfg$scaled_max_outer)) {
      # sparseVAR uses RSS/N; scaled-Lasso subproblem uses 2*lambda0*sigma.
      lam_book <- 2 * lambda0 * sj
      bj <- fista_backend(
        Xs, matrix(Yc[, j], ncol = 1), w,
        lam_book, NULL, cfg
      )
      beta <- as.numeric(bj) / xsd
      r <- Yc[, j] - Xc %*% beta
      sj_new <- sqrt(mean(r^2))
      if (!is.finite(sj_new) || sj_new < 1e-10) sj_new <- 1e-10
      rel <- abs(sj_new - sj) / max(sj, 1e-10)
      sj <- sj_new
      if (rel < cfg$scaled_tol_outer) break
    }
    B[, j] <- beta
    sigma[j] <- sj
    nit[j] <- it
  }

  list(Phi_hat = t(B), lambda0 = lambda0, sigma = sigma, iterations = nit)
}

fit_method_family_vr <- function(X, Y, tuning, cfg) {
  t0 <- proc.time()[3L]
  base <- fit_vr_fista(X, Y, tuning, cfg)

  L <- base$Phi_hat_lasso
  A <- base$Phi_hat_Alasso
  Lref <- vr_refit(X, Y, L)
  Aref <- vr_refit(X, Y, A)

  Lth <- vr_threshold(
    L,
    eta = cfg$threshold_factor * base$lambda_lasso,
    cfg = cfg
  )
  Ath <- vr_threshold(
    A,
    eta = cfg$threshold_factor * base$lambda_adalasso,
    cfg = cfg
  )

  ScObj <- vr_scaled_lasso(X, Y, cfg)
  runtime <- proc.time()[3L] - t0

  list(
    methods = list(
      L = L,
      `L-Th` = Lth,
      `L-Ref` = Lref,
      `L-Rlx` = cfg$relax_gamma * L + (1 - cfg$relax_gamma) * Lref,
      `L-Des` = vr_debias(X, Y, L),
      A = A,
      `A-Th` = Ath,
      `A-Ref` = Aref,
      `A-Rlx` = cfg$relax_gamma * A + (1 - cfg$relax_gamma) * Aref,
      `A-Des` = vr_debias(X, Y, A),
      ScL = ScObj$Phi_hat
    ),
    base = base,
    scaled = ScObj,
    runtime = runtime
  )
}

forecast_rmfe_vr <- function(Phi, X_test, Y_test) {
  pred <- X_test %*% t(Phi)
  sqrt(mean(rowSums((Y_test - pred)^2)))
}

binary_support_metrics <- function(S0, Sh) {
  TP <- sum(Sh & S0)
  FP <- sum(Sh & !S0)
  FN <- sum(!Sh & S0)
  TN <- sum(!Sh & !S0)
  precision <- if (TP + FP > 0) TP / (TP + FP) else 0
  recall <- if (TP + FN > 0) TP / (TP + FN) else 0
  den <- 2 * TP + FP + FN
  F1 <- if (den > 0) 2 * TP / den else NA_real_
  list(
    TP = TP, FP = FP, FN = FN, TN = TN,
    precision = precision, recall = recall, F1 = F1
  )
}

metrics_vr <- function(Phi, Atrue, X_test, Y_test, cfg) {
  d <- nrow(Atrue)
  if (is.null(Phi) || !all(dim(Phi) == dim(Atrue)) || any(!is.finite(Phi))) {
    stop("Invalid/non-finite VR coefficient matrix.")
  }

  S0 <- abs(Atrue) > cfg$support_tol
  Sh <- abs(Phi) > cfg$support_tol
  sm <- binary_support_metrics(S0, Sh)

  off <- row(S0) != col(S0)
  S0_off <- S0[off]
  Sh_off <- Sh[off]
  sm_off <- binary_support_metrics(S0_off, Sh_off)
  # If the true off-diagonal support is empty, off-F1 is not an informative
  # recall/precision score.  Keep NA instead of forcing it to zero.
  if (sum(S0_off) == 0L) sm_off$F1 <- NA_real_

  data.frame(
    frob = sqrt(sum((Phi - Atrue)^2)),
    rmfe = forecast_rmfe_vr(Phi, X_test, Y_test),
    hamming = (sm$FP + sm$FN) / (d * d),
    F1 = if (is.na(sm$F1)) 0 else sm$F1,
    precision = sm$precision,
    recall = sm$recall,
    TP = sm$TP,
    FP = sm$FP,
    FN = sm$FN,
    TN = sm$TN,
    support_true = sum(S0),
    support_hat = sum(Sh),
    off_hamming = (sm_off$FP + sm_off$FN) / sum(off),
    off_F1 = sm_off$F1,
    off_precision = sm_off$precision,
    off_recall = sm_off$recall,
    off_TP = sm_off$TP,
    off_FP = sm_off$FP,
    off_FN = sm_off$FN,
    off_TN = sm_off$TN,
    off_support_true = sum(S0_off),
    off_support_hat = sum(Sh_off),
    stringsAsFactors = FALSE
  )
}

require_vr_packages <- function(cfg) {
  required <- c("sparseVAR", "MASS", "igraph", "foreach", "doParallel")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "))

  got <- as.character(utils::packageVersion("sparseVAR"))
  if (isTRUE(cfg$strict_package_version) && got != cfg$required_sparseVAR) {
    stop("Required sparseVAR ", cfg$required_sparseVAR, "; installed ", got)
  }
  ns <- asNamespace("sparseVAR")
  if (!exists("fista_lasso_multi_cpp", envir = ns, inherits = FALSE)) {
    stop("sparseVAR FISTA backend fista_lasso_multi_cpp is missing.")
  }
  invisible(TRUE)
}

flatten_config <- function(cfg) {
  do.call(rbind, lapply(names(cfg), function(nm) {
    v <- cfg[[nm]]
    val <- if (is.null(v)) "NULL" else paste(as.character(v), collapse = "|")
    data.frame(key = nm, value = val, stringsAsFactors = FALSE)
  }))
}

package_manifest <- function() {
  pkgs <- c("sparseVAR", "MASS", "igraph", "foreach", "doParallel")
  vers <- vapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
  }, character(1))
  data.frame(
    item = c("R", "platform", pkgs),
    version = c(as.character(getRversion()), R.version$platform, unname(vers)),
    stringsAsFactors = FALSE
  )
}

code_manifest <- function(project_dir) {
  files <- sort(list.files(project_dir, pattern = "\\.R$", full.names = TRUE))
  data.frame(
    file = basename(files),
    md5 = unname(tools::md5sum(files)),
    stringsAsFactors = FALSE
  )
}
