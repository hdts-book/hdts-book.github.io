# Same DGPs and seed formulas as the final Chapter 1 VAR methodology simulation.

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
    if (d >= 2L) for (i in seq_len(d - 1L)) {
      M[i, i + 1L] <- stats::runif(1L, cfg$dgp2_b1, cfg$dgp2_b2)
      M[i + 1L, i] <- stats::runif(1L, cfg$dgp2_b1, cfg$dgp2_b2)
    }
  } else if (dgp == "dgp3") {
    bs <- make_block_sizes(d, cfg$dgp3_K)
    P <- matrix(cfg$dgp3_p_out, cfg$dgp3_K, cfg$dgp3_K)
    diag(P) <- cfg$dgp3_p_in
    g <- igraph::sample_sbm(n = d, pref.matrix = P, block.sizes = bs,
                            directed = TRUE, loops = FALSE)
    M <- t(as.matrix(igraph::as_adj(g, sparse = FALSE)))
    diag(M) <- 1
  } else stop("Unknown DGP: ", dgp)
  sm <- largest_singular_value(M)
  if (!is.finite(sm) || sm <= 0) stop("Degenerate DGP matrix.")
  cfg$alpha * M / sm
}

simulate_var1 <- function(A, T_sample, rep_seed, cfg) {
  set_rng_seed(rep_seed, cfg)
  d <- nrow(A)
  TT <- cfg$burn + T_sample + cfg$T_test
  eps <- t(MASS::mvrnorm(TT, mu = rep(0, d), Sigma = diag(d)))
  Y <- matrix(0, d, TT)
  for (tt in 2:TT) Y[, tt] <- A %*% Y[, tt - 1L] + eps[, tt]
  Y <- Y[, (cfg$burn + 1L):TT, drop = FALSE]
  list(
    sample = Y[, seq_len(T_sample), drop = FALSE],
    test = Y[, T_sample + seq_len(cfg$T_test), drop = FALSE]
  )
}
