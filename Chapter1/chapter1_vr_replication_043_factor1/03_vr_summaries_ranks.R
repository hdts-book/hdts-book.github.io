# Summaries and ranking tables for VR 0.4.3 outputs.
# Core book rankings remain Frobenius / RMFE / full-matrix Hamming / F1.
# Off-diagonal support metrics are retained as diagnostics because diagonal
# coefficients are structurally active and unpenalized in all three DGPs.

safe_mean <- function(x) {
  if (length(x) == 0L || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (length(x) == 0L || all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (sum(is.finite(x)) < 2L) NA_real_ else stats::sd(x, na.rm = TRUE)
}

summarize_vr_cells <- function(raw, cfg, project_dir, method_order) {
  out_root <- file.path(project_dir, cfg$output_dir)

  vars <- c(
    "frob", "rmfe", "hamming", "F1",
    "precision", "recall", "TP", "FP", "FN", "TN",
    "support_true", "support_hat",
    "off_hamming", "off_F1", "off_precision", "off_recall",
    "off_TP", "off_FP", "off_FN", "off_TN",
    "off_support_true", "off_support_hat",
    "lambda_lasso", "lambda_adaptive", "scaled_lambda0",
    "base_lasso_support", "base_adaptive_support", "runtime_fit_seconds"
  )

  key <- interaction(
    raw$snr, raw$tuning, raw$dgp, raw$d, raw$T_sample, raw$method,
    drop = TRUE, lex.order = TRUE
  )

  parts <- lapply(split(raw, key), function(z) {
    r <- z[1, c("snr", "tuning", "dgp", "d", "T_sample", "method"), drop = FALSE]
    r$nrep_eff <- nrow(z)
    for (v in vars) {
      x <- z[[v]]
      r[[paste0(v, "_mean")]] <- safe_mean(x)
      r[[paste0(v, "_median")]] <- safe_median(x)
      r[[paste0(v, "_sd")]] <- safe_sd(x)
    }
    r
  })

  cell <- do.call(rbind, parts)
  rownames(cell) <- NULL
  cell <- cell[order(
    match(cell$snr, cfg$snr_grid),
    match(cell$tuning, cfg$tunings),
    match(cell$dgp, cfg$dgps),
    cell$d, cell$T_sample,
    match(cell$method, method_order)
  ), , drop = FALSE]

  expected <- length(cfg$snr_grid) * length(cfg$tunings) *
    length(cfg$dgps) * length(cfg$d_grid) * length(cfg$T_grid) * length(method_order)
  if (nrow(cell) != expected) stop("Cell-summary row-count mismatch.")
  if (any(cell$nrep_eff != cfg$nrep)) stop("Incomplete VR cell summary.")

  utils::write.csv(cell, file.path(out_root, "summary", "cell_performance.csv"), row.names = FALSE)
  cell
}

rank_or_na <- function(x, decreasing = FALSE) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rank(if (decreasing) -x else x, ties.method = "average", na.last = "keep")
}

make_vr_cell_ranks <- function(cell, cfg, project_dir, method_order) {
  out_root <- file.path(project_dir, cfg$output_dir)
  key <- interaction(
    cell$snr, cell$tuning, cell$dgp, cell$d, cell$T_sample,
    drop = TRUE, lex.order = TRUE
  )

  parts <- lapply(split(cell, key), function(z) {
    if (!setequal(z$method, method_order)) stop("Incomplete method set while ranking VR cell.")
    z <- z[match(method_order, z$method), , drop = FALSE]
    data.frame(
      snr = z$snr,
      tuning = z$tuning,
      dgp = z$dgp,
      d = z$d,
      T_sample = z$T_sample,
      method = z$method,
      rank_frob = rank_or_na(z$frob_median),
      rank_rmfe = rank_or_na(z$rmfe_median),
      rank_hamming = rank_or_na(z$hamming_median),
      rank_F1 = rank_or_na(z$F1_median, decreasing = TRUE),
      rank_off_hamming = rank_or_na(z$off_hamming_median),
      rank_off_F1 = rank_or_na(z$off_F1_median, decreasing = TRUE),
      stringsAsFactors = FALSE
    )
  })

  rr <- do.call(rbind, parts)
  rownames(rr) <- NULL
  rr <- rr[order(
    match(rr$snr, cfg$snr_grid),
    match(rr$tuning, cfg$tunings),
    match(rr$dgp, cfg$dgps),
    rr$d, rr$T_sample,
    match(rr$method, method_order)
  ), , drop = FALSE]

  expected <- length(cfg$snr_grid) * length(cfg$tunings) * 108L * length(method_order)
  if (nrow(rr) != expected) stop("Cell-rank row-count mismatch.")

  utils::write.csv(rr, file.path(out_root, "ranks", "cell_ranks.csv"), row.names = FALSE)
  utils::write.csv(rr[rr$snr == 1, , drop = FALSE], file.path(out_root, "ranks", "cell_ranks_snr1.csv"), row.names = FALSE)
  utils::write.csv(rr[rr$snr == 2, , drop = FALSE], file.path(out_root, "ranks", "cell_ranks_snr2.csv"), row.names = FALSE)
  rr
}

aggregate_vr_ranks <- function(rr, cfg, project_dir, method_order) {
  out_root <- file.path(project_dir, cfg$output_dir)

  aggregate_one <- function(group_cols, filename) {
    split_cols <- c(group_cols, "method")
    key <- interaction(rr[split_cols], drop = TRUE, lex.order = TRUE)
    pieces <- lapply(split(rr, key), function(z) {
      r <- z[1, split_cols, drop = FALSE]
      r$n_cells <- nrow(z)
      r$avg_rank_frob <- safe_mean(z$rank_frob)
      r$avg_rank_rmfe <- safe_mean(z$rank_rmfe)
      r$avg_rank_hamming <- safe_mean(z$rank_hamming)
      r$avg_rank_F1 <- safe_mean(z$rank_F1)
      r$avg_rank_off_hamming <- safe_mean(z$rank_off_hamming)
      r$avg_rank_off_F1 <- safe_mean(z$rank_off_F1)
      r
    })
    ans <- do.call(rbind, pieces)
    ans <- ans[order(
      match(ans$snr, cfg$snr_grid),
      match(ans$tuning, cfg$tunings),
      if ("dgp" %in% names(ans)) match(ans$dgp, cfg$dgps) else rep(0L, nrow(ans)),
      if ("d" %in% names(ans)) ans$d else rep(0L, nrow(ans)),
      if ("T_sample" %in% names(ans)) ans$T_sample else rep(0L, nrow(ans)),
      match(ans$method, method_order)
    ), , drop = FALSE]
    rownames(ans) <- NULL
    utils::write.csv(ans, file.path(out_root, "ranks", filename), row.names = FALSE)
    ans
  }

  overall <- aggregate_one(c("snr", "tuning"), "average_ranks_overall.csv")
  aggregate_one(c("snr", "tuning", "dgp"), "average_ranks_by_dgp.csv")
  aggregate_one(c("snr", "tuning", "d"), "average_ranks_by_d.csv")
  aggregate_one(c("snr", "tuning", "T_sample"), "average_ranks_by_T.csv")
  overall
}

aggregate_vr_performance <- function(cell, cfg, project_dir, method_order) {
  out_root <- file.path(project_dir, cfg$output_dir)
  perf_cols <- c("frob", "rmfe", "hamming", "F1", "off_hamming", "off_F1")

  aggregate_one <- function(group_cols, filename) {
    split_cols <- c(group_cols, "method")
    key <- interaction(cell[split_cols], drop = TRUE, lex.order = TRUE)
    pieces <- lapply(split(cell, key), function(z) {
      r <- z[1, split_cols, drop = FALSE]
      r$n_cells <- nrow(z)
      for (v in perf_cols) {
        r[[paste0("avg_cell_mean_", v)]] <- safe_mean(z[[paste0(v, "_mean")]])
        r[[paste0("avg_cell_median_", v)]] <- safe_mean(z[[paste0(v, "_median")]])
      }
      r
    })
    ans <- do.call(rbind, pieces)
    ans <- ans[order(
      match(ans$snr, cfg$snr_grid),
      match(ans$tuning, cfg$tunings),
      if ("dgp" %in% names(ans)) match(ans$dgp, cfg$dgps) else rep(0L, nrow(ans)),
      if ("d" %in% names(ans)) ans$d else rep(0L, nrow(ans)),
      if ("T_sample" %in% names(ans)) ans$T_sample else rep(0L, nrow(ans)),
      match(ans$method, method_order)
    ), , drop = FALSE]
    rownames(ans) <- NULL
    utils::write.csv(ans, file.path(out_root, "summary", filename), row.names = FALSE)
    ans
  }

  overall <- aggregate_one(c("snr", "tuning"), "average_performance_overall.csv")
  aggregate_one(c("snr", "tuning", "dgp"), "average_performance_by_dgp.csv")
  aggregate_one(c("snr", "tuning", "d"), "average_performance_by_d.csv")
  aggregate_one(c("snr", "tuning", "T_sample"), "average_performance_by_T.csv")
  overall
}

make_vr_rank_table_wide <- function(overall, method_order) {
  get_block <- function(snr, tuning) {
    z <- overall[overall$snr == snr & overall$tuning == tuning, , drop = FALSE]
    z[match(method_order, z$method), , drop = FALSE]
  }
  blocks <- list(
    VR_Th = get_block(1, "theory"),
    VR_CV = get_block(1, "cv"),
    VR_SNR2_Th = get_block(2, "theory"),
    VR_SNR2_CV = get_block(2, "cv")
  )
  out <- data.frame(method = method_order, stringsAsFactors = FALSE)
  for (nm in names(blocks)) {
    z <- blocks[[nm]]
    for (metric in c("frob", "rmfe", "hamming", "F1")) {
      out[[paste0(nm, "_rank_", metric)]] <- z[[paste0("avg_rank_", metric)]]
    }
  }
  out
}

write_vr_wide_tables <- function(overall_rank, cfg, project_dir, method_order) {
  out_root <- file.path(project_dir, cfg$output_dir)
  wide <- make_vr_rank_table_wide(overall_rank, method_order)
  utils::write.csv(wide, file.path(out_root, "ranks", "rank_table_overall_wide.csv"), row.names = FALSE)
  invisible(wide)
}
