# Historical Chapter 1 DGP1--DGP3 and deterministic VAR(1) simulation.
# Definitions and seed formulas are intentionally unchanged from 0.4.2.

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
      n = d, pref.matrix = P, block.sizes = bs,
      directed = TRUE, loops = FALSE
    )
    M <- t(as.matrix(igraph::as_adj(g, sparse = FALSE)))
    diag(M) <- 1
  } else {
    stop("Unknown DGP: ", dgp)
  }

  sm <- largest_singular_value(M)
  if (!is.finite(sm) || sm <= 0) stop("Degenerate DGP matrix.")
  cfg$alpha * M / sm
}

simulate_var1 <- function(A, T_sample, rep_seed, cfg) {
  set_rng_seed(rep_seed, cfg)
  d <- nrow(A)
  TT <- cfg$burn + T_sample + cfg$T_test
  # Innovation covariance is exactly I_d, matching the book-theory naive loss.
  eps <- t(MASS::mvrnorm(TT, mu = rep(0, d), Sigma = diag(d)))
  Y <- matrix(0, d, TT)
  for (tt in 2:TT) Y[, tt] <- A %*% Y[, tt - 1L] + eps[, tt]
  Y <- Y[, (cfg$burn + 1L):TT, drop = FALSE]
  list(
    sample = Y[, seq_len(T_sample), drop = FALSE],
    test = Y[, T_sample + seq_len(cfg$T_test), drop = FALSE]
  )
}

truth_matrix_long <- function(cfg) {
  pieces <- list(); kk <- 1L
  for (dgp in cfg$dgps) {
    for (idx_d in seq_along(cfg$d_grid)) {
      d <- cfg$d_grid[idx_d]
      A <- make_A_true(dgp, d, idx_d, cfg)
      if (any(abs(diag(A)) <= cfg$support_tol)) {
        stop("Structural diagonal audit failed for ", dgp, ", d=", d)
      }
      zz <- expand.grid(row = seq_len(d), col = seq_len(d), KEEP.OUT.ATTRS = FALSE)
      zz$dgp <- dgp; zz$idx_d <- idx_d; zz$d <- d
      zz$truth_seed <- truth_seed(dgp, idx_d, cfg$truth_seed_base)
      zz$value <- A[cbind(zz$row, zz$col)]
      zz$rho_true <- spectral_radius(A)
      zz$support_true <- sum(abs(A) > cfg$support_tol)
      off <- row(A) != col(A)
      zz$support_true_offdiag <- sum(abs(A[off]) > cfg$support_tol)
      pieces[[kk]] <- zz[, c(
        "dgp", "idx_d", "d", "truth_seed", "row", "col", "value",
        "rho_true", "support_true", "support_true_offdiag"
      )]
      kk <- kk + 1L
    }
  }
  do.call(rbind, pieces)
}
