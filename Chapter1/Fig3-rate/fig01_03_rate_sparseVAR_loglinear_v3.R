# =============================================================================
# Figure 1.3 — Sparse VAR convergence rate
# sparseVAR package version
#
# Basis:
#   sparseVAR 0.3.0
#
# This script:
#   * depends on sparseVAR, not HDTS;
#   * simulates the random-graph stable VAR(1) used in Chapter 1;
#   * fits a plain LASSO using sparseVAR's compiled FISTA backend;
#   * uses theoretical tuning lambda_N = c * sqrt(log(q)/N);
#   * converts lambda_N to the native sparseVAR FISTA scale correctly;
#   * computes Monte Carlo mean squared l2/Frobenius error directly;
#   * plots
#
#       left:  E ||beta_hat - beta*||_2^2  versus T
#
#       right: log{N/(s log q)}
#              versus
#              log ||beta_hat - beta*||_2^2
#
#     with a theoretical slope -1 reference line.
#
# IMPORTANT sparseVAR scaling
# ---------------------------
# sparseVAR's C++ FISTA backend solves
#
#   0.5 ||Y - X B||_F^2
#       + lambda_native * sum_ij w_ij |B_ij|.
#
# The Chapter-1 normalized LASSO criterion is
#
#   (1/(2N)) ||Y - X B||_F^2
#       + lambda_N ||B||_1.
#
# Multiplying the latter by N gives the sparseVAR native scale:
#
#   lambda_native = N * lambda_N.
#
# For lambda_N = c * sqrt(log(q)/N),
#
#   lambda_native = c * sqrt(N log q).
#
# =============================================================================


# =============================================================================
# Source-safe entry point
# =============================================================================
# All executable code is inside main(). This is deliberate: calling on.exit()
# at the top level of a script that is executed with source() can interfere
# with source()'s own connection cleanup on Windows.
main <- function() {


  # -----------------------------------------------------------------------------
  # Settings
  # -----------------------------------------------------------------------------
  ndim <- c(10, 15, 20, 25, 30)
  nT   <- c(100, 200, 300, 500, 1000, 2000)

  p <- 1L

  # Random-graph DGP
  prop   <- 0.10
  c1     <- 0.90
  nalpha <- 0.90

  # Monte Carlo
  nrep <- 100L
  # Windows/PSOCK: a moderate default is more robust than opening 15 workers.
  # Increase this manually if desired after confirming the script runs.
  ncore <- min(
    8L,
    max(1L, parallel::detectCores(logical = FALSE) - 1L)
  )

  base_seed <- 20260907L

  # Theoretical penalty multiplier:
  # lambda_N = lambda_constant * sqrt(log(q)/N)
  lambda_constant <- 1.0

  # Retain the convention used in the existing Figure 1.3 code:
  # TRUE = own-lag diagonal coefficients are not penalized.
  diagTF <- TRUE

  # FISTA controls
  max_iter <- 3000L
  tol <- 1e-7

  # Save replication-level errors?
  # TRUE is useful for diagnostics but makes the RDS larger.
  save_replication_errors <- TRUE

  # Output
  out_rds   <- "fig01_03_sparseVAR_results.rds"
  out_rdata <- "fig01_03_sparseVAR_results.RData"

  fig_pdf <- "Ch1-sVAR-error-rates-sparseVAR.pdf"
  fig_png <- "Ch1-sVAR-error-rates-sparseVAR.png"
  slope_csv <- "Ch1-sVAR-error-rates-sparseVAR-slopes.csv"

  # -----------------------------------------------------------------------------
  # Package
  # -----------------------------------------------------------------------------
  if (!requireNamespace("sparseVAR", quietly = TRUE)) {
    stop(
      "Package 'sparseVAR' is required.\n",
      "Install sparseVAR_0.3.0.tar.gz first, for example:\n",
      "  install.packages('sparseVAR_0.3.0.tar.gz', repos = NULL, type = 'source')"
    )
  }

  if (packageVersion("sparseVAR") < package_version("0.3.0")) {
    stop("This script requires sparseVAR >= 0.3.0.")
  }

  required_plot_packages <- c(
    "ggplot2",
    "RColorBrewer",
    "ggpubr",
    "foreach",
    "doParallel"
  )

  missing_plot_packages <- required_plot_packages[
    !vapply(
      required_plot_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(missing_plot_packages)) {
    stop(
      "Install required package(s): ",
      paste(missing_plot_packages, collapse = ", ")
    )
  }

  suppressPackageStartupMessages({
    library(sparseVAR)
    library(ggplot2)
    library(RColorBrewer)
    library(ggpubr)
    library(foreach)
    library(doParallel)
  })

  # -----------------------------------------------------------------------------
  # Check the package backends required for fixed-theory FISTA.
  #
  # sparseVAR 0.3.0 exports sVAR_adalasso_fista(), but that wrapper chooses
  # lambda by blocked CV and does not expose a fixed theoretical lambda.
  # For this convergence-rate experiment we therefore call the package's
  # compiled FISTA backend directly through its namespace.
  # -----------------------------------------------------------------------------
  svar_ns <- asNamespace("sparseVAR")

  needed_internal <- c(
    "VAR_sim_cpp",
    "create_var_design_matrix_cpp",
    "fista_lasso_multi_cpp"
  )

  missing_internal <- needed_internal[
    !vapply(
      needed_internal,
      exists,
      logical(1),
      envir = svar_ns,
      inherits = FALSE
    )
  ]

  if (length(missing_internal)) {
    stop(
      "The installed sparseVAR package is missing required compiled backends: ",
      paste(missing_internal, collapse = ", "),
      ". Reinstall/rebuild sparseVAR 0.3.0."
    )
  }

  VAR_sim_cpp <- get(
    "VAR_sim_cpp",
    envir = svar_ns,
    inherits = FALSE
  )

  create_var_design_matrix_cpp <- get(
    "create_var_design_matrix_cpp",
    envir = svar_ns,
    inherits = FALSE
  )

  fista_lasso_multi_cpp <- get(
    "fista_lasso_multi_cpp",
    envir = svar_ns,
    inherits = FALSE
  )

  # -----------------------------------------------------------------------------
  # Random-graph stable VAR(1) coefficient matrix
  #
  # Generate ONE Phi* per dimension and reuse it across all T.
  # This is preferable for a convergence-rate curve: as T changes, the DGP is
  # held fixed and only sample size changes.
  # -----------------------------------------------------------------------------
  make_random_graph_var1 <- function(
      d,
      prop = 0.10,
      c1 = 0.90,
      nalpha = 0.90) {

    A <- matrix(
      sample(
        c(0, 1),
        d * d,
        replace = TRUE,
        prob = c(1 - prop, prop)
      ),
      nrow = d,
      ncol = d
    )

    diag(A) <- 0

    # Symmetric adjacency pattern
    A <- 1L * ((A + t(A)) > 0)

    # Include own-lag effects
    A2 <- A + diag(d)

    nk <- rowSums(A2 != 0)

    B1 <- c1 * diag(1 / nk, d, d) %*% A2

    # Chapter-1 scaling:
    # ||Phi*||_2 = nalpha < 1, hence rho(Phi*) < 1.
    sigma_max <- max(svd(B1, nu = 0, nv = 0)$d)

    Phi <- (nalpha / sigma_max) * B1

    list(
      Phi = Phi,
      support_size = sum(Phi != 0),
      spectral_norm = max(svd(Phi, nu = 0, nv = 0)$d),
      spectral_radius = max(Mod(eigen(Phi, only.values = TRUE)$values))
    )
  }

  # -----------------------------------------------------------------------------
  # Fixed-lambda sparseVAR FISTA LASSO
  #
  # Input Yt: d x T
  # Output Phi_hat: d x d for VAR(1)
  # -----------------------------------------------------------------------------
  fit_sparseVAR_fista_theory <- function(
      Yt,
      p = 1L,
      lambda_N,
      diagTF = TRUE,
      max_iter = 3000L,
      tol = 1e-7) {

    Yt <- as.matrix(Yt)

    k  <- nrow(Yt)
    TT <- ncol(Yt)

    if (p < 1L || p >= TT) {
      stop("Require 1 <= p < T.")
    }

    # Package wrappers use no intercept and center each series first.
    Yt_centered <- Yt - rowMeans(Yt)

    # sparseVAR C++ design helper takes observations in rows.
    des <- create_var_design_matrix_cpp(
      X = t(Yt_centered),
      p = as.integer(p)
    )

    X <- as.matrix(des$X_design)  # N x (k p)
    Y <- as.matrix(des$Y_resp)    # N x k

    N <- nrow(X)

    # Weights have dimension (k p) x k.
    W <- matrix(
      1,
      nrow = k * p,
      ncol = k
    )

    if (isTRUE(diagTF)) {
      for (ell in seq_len(p)) {
        for (j in seq_len(k)) {
          row_id <- (ell - 1L) * k + j
          W[row_id, j] <- 0
        }
      }
    }

    # -----------------------------------------------------------
    # Critical conversion:
    #
    # book:
    #   (1/(2N)) RSS + lambda_N penalty
    #
    # sparseVAR FISTA:
    #   0.5 RSS + lambda_native penalty
    #
    # so lambda_native = N * lambda_N.
    # -----------------------------------------------------------
    lambda_native <- N * lambda_N

    B_hat <- fista_lasso_multi_cpp(
      X = X,
      Y = Y,
      weights = W,
      lambda = lambda_native,
      SigmaInv = matrix(0, 0, 0),
      max_iter = as.integer(max_iter),
      tol = tol
    )

    # sparseVAR design is Y = X B + U.
    # Chapter VAR convention is X_t = Phi X_{t-1} + eps_t,
    # therefore Phi_hat = B_hat'.
    Phi_hat <- t(B_hat)

    list(
      Phi_hat = Phi_hat,
      lambda_N = lambda_N,
      lambda_native = lambda_native,
      N = N
    )
  }

  # -----------------------------------------------------------------------------
  # Generate fixed DGPs by dimension
  # -----------------------------------------------------------------------------
  set.seed(base_seed)

  dgp <- vector("list", length(ndim))
  names(dgp) <- as.character(ndim)

  for (ii in seq_along(ndim)) {

    d <- ndim[ii]

    set.seed(base_seed + 1000L * ii)

    dgp[[ii]] <- make_random_graph_var1(
      d = d,
      prop = prop,
      c1 = c1,
      nalpha = nalpha
    )

    cat(
      sprintf(
        "[DGP] d=%d, s=%d, ||Phi||2=%.4f, rho(Phi)=%.4f\n",
        d,
        dgp[[ii]]$support_size,
        dgp[[ii]]$spectral_norm,
        dgp[[ii]]$spectral_radius
      )
    )
  }

  # -----------------------------------------------------------------------------
  # Parallel backend
  # -----------------------------------------------------------------------------
  cl <- parallel::makeCluster(
    ncore,
    type = "PSOCK",
    outfile = ""
  )

  doParallel::registerDoParallel(cl)

  # IMPORTANT:
  # This on.exit() is now inside main(), so it cannot interfere with
  # source()'s own file connection.
  on.exit(
    try(parallel::stopCluster(cl), silent = TRUE),
    add = TRUE
  )

  # Confirm that every worker can load sparseVAR before the simulation starts.
  worker_check <- parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages(
      library(sparseVAR)
    )
    list(
      pid = Sys.getpid(),
      sparseVAR = as.character(packageVersion("sparseVAR"))
    )
  })

  cat(
    sprintf(
      "[cluster] %d PSOCK workers started; sparseVAR versions: %s\n",
      length(worker_check),
      paste(
        unique(vapply(worker_check, `[[`, character(1), "sparseVAR")),
        collapse = ", "
      )
    )
  )

  # -----------------------------------------------------------------------------
  # Storage
  # -----------------------------------------------------------------------------
  nd <- length(ndim)
  nt <- length(nT)

  mse <- matrix(
    NA_real_,
    nrow = nd,
    ncol = nt,
    dimnames = list(
      paste0("d=", ndim),
      paste0("T=", nT)
    )
  )

  n_success <- matrix(
    0L,
    nrow = nd,
    ncol = nt,
    dimnames = dimnames(mse)
  )

  lambda_book <- matrix(
    NA_real_,
    nrow = nd,
    ncol = nt,
    dimnames = dimnames(mse)
  )

  lambda_native <- matrix(
    NA_real_,
    nrow = nd,
    ncol = nt,
    dimnames = dimnames(mse)
  )

  replication_errors <- if (save_replication_errors) {
    vector("list", nd * nt)
  } else {
    NULL
  }

  if (save_replication_errors) {
    dim(replication_errors) <- c(nd, nt)
  }

  # -----------------------------------------------------------------------------
  # Monte Carlo
  # -----------------------------------------------------------------------------
  cell_id <- 0L

  for (ii in seq_along(ndim)) {

    d <- ndim[ii]

    Phi_true <- dgp[[ii]]$Phi
    s_true   <- dgp[[ii]]$support_size

    Sigma_e <- diag(d)

    q <- p * d^2

    for (jj in seq_along(nT)) {

      TT <- nT[jj]
      N  <- TT - p

      # Chapter theoretical rate
      lam_N <- lambda_constant * sqrt(
        log(q) / N
      )

      # sparseVAR native unnormalized objective
      lam_native <- N * lam_N

      cell_id <- cell_id + 1L

      # Reproducible PSOCK streams by (d,T) cell
      parallel::clusterSetRNGStream(
        cl,
        iseed = base_seed + 100000L + cell_id
      )

      sqerr <- foreach(
        r = seq_len(nrep),
        .packages = "sparseVAR",
        .errorhandling = "remove"
      ) %dopar% {

        ns <- asNamespace("sparseVAR")

        sim_fun <- get(
          "VAR_sim_cpp",
          envir = ns,
          inherits = FALSE
        )

        design_fun <- get(
          "create_var_design_matrix_cpp",
          envir = ns,
          inherits = FALSE
        )

        fista_fun <- get(
          "fista_lasso_multi_cpp",
          envir = ns,
          inherits = FALSE
        )

        # Simulate using sparseVAR's compiled VAR simulator.
        Yt <- sim_fun(
          T = as.integer(TT),
          A = Phi_true,
          Sigma = Sigma_e,
          burn = 500L
        )

        # Center exactly as sparseVAR's high-level FISTA wrapper does.
        Ytc <- Yt - rowMeans(Yt)

        des <- design_fun(
          X = t(Ytc),
          p = as.integer(p)
        )

        X <- as.matrix(des$X_design)
        Y <- as.matrix(des$Y_resp)

        W <- matrix(
          1,
          nrow = d * p,
          ncol = d
        )

        if (isTRUE(diagTF)) {
          for (ell in seq_len(p)) {
            for (j in seq_len(d)) {
              W[(ell - 1L) * d + j, j] <- 0
            }
          }
        }

        Bhat <- fista_fun(
          X = X,
          Y = Y,
          weights = W,
          lambda = lam_native,
          SigmaInv = matrix(0, 0, 0),
          max_iter = as.integer(max_iter),
          tol = tol
        )

        Phi_hat <- t(Bhat)

        # Squared vectorized l2 error = squared Frobenius error.
        sum(
          (Phi_hat - Phi_true)^2
        )
      }

      sqerr <- unlist(
        sqerr,
        use.names = FALSE
      )

      if (!length(sqerr)) {
        stop(
          sprintf(
            "All replications failed for d=%d, T=%d.",
            d,
            TT
          )
        )
      }

      mse[ii, jj] <- mean(
        sqerr,
        na.rm = TRUE
      )

      n_success[ii, jj] <- length(sqerr)

      lambda_book[ii, jj] <- lam_N
      lambda_native[ii, jj] <- lam_native

      if (save_replication_errors) {
        replication_errors[[ii, jj]] <- sqerr
      }

      cat(
        sprintf(
          paste0(
            "[done] d=%d, T=%d, N=%d, s=%d, ",
            "lambda_N=%.5f, lambda_native=%.3f, ",
            "MSE=%.5f, success=%d/%d\n"
          ),
          d,
          TT,
          N,
          s_true,
          lam_N,
          lam_native,
          mse[ii, jj],
          n_success[ii, jj],
          nrep
        )
      )
    }
  }

  parallel::stopCluster(cl)

  # Remove only this main() cleanup handler after the explicit successful stop.
  on.exit(NULL, add = FALSE)
  cl <- NULL

  # -----------------------------------------------------------------------------
  # Save simulation results
  # -----------------------------------------------------------------------------
  out <- list(
    package = "sparseVAR",
    package_version = as.character(
      packageVersion("sparseVAR")
    ),
    ndim = ndim,
    nT = nT,
    p = p,
    prop = prop,
    c1 = c1,
    nalpha = nalpha,
    nrep = nrep,
    diagTF = diagTF,
    lambda_constant = lambda_constant,
    dgp = dgp,
    mse = mse,
    n_success = n_success,
    lambda_book = lambda_book,
    lambda_native = lambda_native,
    replication_errors = replication_errors
  )

  saveRDS(
    out,
    file = out_rds
  )

  save(
    out,
    file = out_rdata
  )

  # -----------------------------------------------------------------------------
  # Plot data
  # -----------------------------------------------------------------------------
  ff <- expand.grid(
    T = nT,
    d = ndim,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  # t(mse): rows T, columns d.
  ff$Value <- as.vector(
    t(mse)
  )

  ff$Series <- factor(
    sprintf("$d=%d$", ff$d),
    levels = sprintf(
      "$d=%d$",
      ndim
    )
  )

  d_index <- match(
    ff$d,
    ndim
  )

  ff$s <- vapply(
    d_index,
    function(i) dgp[[i]]$support_size,
    numeric(1)
  )

  ff$N <- ff$T - p
  ff$q <- p * ff$d^2

  ff$Nscaled <- ff$N / (
    ff$s * log(ff$q)
  )

  if (any(
    !is.finite(ff$Value) |
      ff$Value <= 0
  )) {
    stop(
      "Positive MSE values are required for log(error)."
    )
  }

  if (any(
    !is.finite(ff$Nscaled) |
      ff$Nscaled <= 0
  )) {
    stop(
      "Positive N/(s log q) values are required."
    )
  }

  # Explicit natural-log transformed variables
  ff$log_Nscaled <- log(
    ff$Nscaled
  )

  ff$log_Error <- log(
    ff$Value
  )

  # -----------------------------------------------------------------------------
  # Palette
  # -----------------------------------------------------------------------------
  coul <- colorRampPalette(
    brewer.pal(
      4,
      "PuOr"
    )
  )(25)

  plot_cols <- coul[
    c(1, 5, 9, 13, 17)
  ]

  names(plot_cols) <- levels(
    ff$Series
  )

  # -----------------------------------------------------------------------------
  # Left panel
  # -----------------------------------------------------------------------------
  p_T <- ggplot(
    ff,
    aes(
      x = T,
      y = Value,
      group = Series,
      color = Series
    )
  ) +
    geom_line(
      linewidth = 0.55
    ) +
    geom_point(
      size = 1.7
    ) +
    scale_color_manual(
      values = plot_cols
    ) +
    labs(
      x = "$T$",
      y = "$\\|\\widehat{\\beta}-\\beta^*\\|_2^2$",
      color = NULL
    ) +
    theme_bw(
      base_size = 10
    ) +
    theme(
      legend.position = "top",
      legend.title = element_blank()
    )

  # -----------------------------------------------------------------------------
  # Right panel
  #
  # log(error) = constant - log{N/(s log q)}
  # theoretical slope = -1
  # -----------------------------------------------------------------------------
  a_ref <- mean(
    ff$log_Error +
      ff$log_Nscaled,
    na.rm = TRUE
  )

  ref_df <- data.frame(
    log_Nscaled = seq(
      min(ff$log_Nscaled),
      max(ff$log_Nscaled),
      length.out = 200
    )
  )

  ref_df$log_Error <- a_ref -
    ref_df$log_Nscaled

  x_lab <- min(
    ff$log_Nscaled
  ) +
    0.68 * diff(
      range(ff$log_Nscaled)
    )

  y_lab <- a_ref -
    x_lab +
    0.10

  p_log <- ggplot(
    ff,
    aes(
      x = log_Nscaled,
      y = log_Error,
      group = Series,
      color = Series
    )
  ) +
    geom_line(
      linewidth = 0.55
    ) +
    geom_point(
      size = 1.7
    ) +
    geom_line(
      data = ref_df,
      aes(
        x = log_Nscaled,
        y = log_Error
      ),
      inherit.aes = FALSE,
      linewidth = 0.65,
      linetype = "dashed",
      color = "black"
    ) +
    annotate(
      "text",
      x = x_lab,
      y = y_lab,
      label = "slope -1",
      hjust = 0,
      size = 3.2
    ) +
    scale_color_manual(
      values = plot_cols
    ) +
    labs(
      x = "$\\log\\{N/(s\\log q)\\}$",
      y = "$\\log\\|\\widehat{\\beta}-\\beta^*\\|_2^2$",
      color = NULL
    ) +
    theme_bw(
      base_size = 10
    ) +
    theme(
      legend.position = "top",
      legend.title = element_blank()
    )

  # -----------------------------------------------------------------------------
  # Descriptive log-linear slopes
  # -----------------------------------------------------------------------------
  slope_by_d <- do.call(
    rbind,
    lapply(
      split(
        ff,
        ff$d
      ),
      function(dd) {

        fit <- lm(
          log_Error ~ log_Nscaled,
          data = dd
        )

        data.frame(
          d = unique(dd$d),
          support_size = unique(dd$s),
          slope = unname(
            coef(fit)[2]
          ),
          r_squared = summary(
            fit
          )$r.squared
        )
      }
    )
  )

  rownames(
    slope_by_d
  ) <- NULL

  fit_all <- lm(
    log_Error ~ log_Nscaled,
    data = ff
  )

  cat(
    "\nFinite-sample slopes by dimension:\n"
  )

  print(
    slope_by_d,
    row.names = FALSE,
    digits = 4
  )

  cat(
    sprintf(
      "\nPooled descriptive slope = %.4f (R^2 = %.4f)\n",
      coef(fit_all)[2],
      summary(fit_all)$r.squared
    )
  )

  cat(
    "Theoretical benchmark: slope = -1.\n\n"
  )

  write.csv(
    slope_by_d,
    file = slope_csv,
    row.names = FALSE
  )

  # -----------------------------------------------------------------------------
  # Combine and save
  # -----------------------------------------------------------------------------
  p_both <- ggarrange(
    p_T,
    p_log,
    ncol = 2,
    nrow = 1,
    common.legend = TRUE,
    legend = "top",
    align = "hv"
  )

  print(
    p_both
  )

  ggsave(
    filename = fig_pdf,
    plot = p_both,
    width = 7.1,
    height = 3.8,
    units = "in"
  )

  ggsave(
    filename = fig_png,
    plot = p_both,
    width = 7.1,
    height = 3.8,
    units = "in",
    dpi = 250
  )

  # Optional TikZ output
  #
  # IMPORTANT: labels above are raw LaTeX strings.  sanitize=FALSE is
  # required so that tikzDevice passes $...$ through unchanged.
  if (
    requireNamespace(
      "tikzDevice",
      quietly = TRUE
    )
  ) {

    tikz_file <-
      "Ch1-sVAR-error-rates-sparseVAR.tex"

    tikzDevice::tikz(
      file = tikz_file,
      width = 7.1,
      height = 3.8,
      standAlone = TRUE,
      sanitize = FALSE
    )

    print(
      p_both
    )

    dev.off()

    message(
      "TikZ source written to: ",
      tikz_file
    )
  }

  cat(
    "\nSaved:\n",
    "  ", out_rds, "\n",
    "  ", out_rdata, "\n",
    "  ", fig_pdf, "\n",
    "  ", fig_png, "\n",
    "  ", slope_csv, "\n",
    sep = ""
  )
}

# Run when sourced or executed via Rscript.
main()
