# Cell ranks and aggregate rank/performance tables.
# Main four book criteria are unchanged.  Off-diagonal support metrics are
# retained as additional diagnostics because all three DGPs have active diagonals.

rank_cell_methods <- function(cell_perf, cfg, method_order) {
  keys <- interaction(cell_perf[c("dgp", "d", "T_sample", "tuning")], drop = TRUE, lex.order = TRUE)
  pieces <- lapply(split(cell_perf, keys), function(a) {
    if (!setequal(a$method, method_order)) stop("Incomplete method set while computing cell ranks.")
    a <- a[match(method_order, a$method), , drop = FALSE]
    a$rank_frob <- rank(a$frob_median, ties.method = "average", na.last = "keep")
    a$rank_rmfe <- rank(a$rmfe_median, ties.method = "average", na.last = "keep")
    a$rank_hamming <- rank(a$hamming_median, ties.method = "average", na.last = "keep")
    a$rank_F1 <- rank(-a$F1_median, ties.method = "average", na.last = "keep")
    a$rank_off_hamming <- rank(a$off_hamming_median, ties.method = "average", na.last = "keep")
    a$rank_off_F1 <- rank(-a$off_F1_median, ties.method = "average", na.last = "keep")
    a
  })
  out <- do.call(rbind, pieces)
  out <- out[
    order(match(out$tuning, cfg$tunings), match(out$dgp, cfg$dgps),
          out$d, out$T_sample, match(out$method, method_order)),
    , drop = FALSE
  ]
  rownames(out) <- NULL
  out
}

aggregate_ranks <- function(cell_ranks, group_cols, method_order) {
  split_cols <- c(group_cols, "method")
  key <- interaction(cell_ranks[split_cols], drop = TRUE, lex.order = TRUE)
  pieces <- lapply(split(cell_ranks, key), function(a) {
    r <- a[1, split_cols, drop = FALSE]
    r$n_cells <- nrow(a)
    for (m in c("frob", "rmfe", "hamming", "F1", "off_hamming", "off_F1")) {
      x <- a[[paste0("rank_", m)]]
      r[[paste0("avg_rank_", m)]] <- if (any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_
      r[[paste0("n_ranked_", m)]] <- sum(is.finite(x))
    }
    r
  })
  out <- do.call(rbind, pieces)
  ord_args <- c(lapply(group_cols, function(g) out[[g]]), list(match(out$method, method_order)))
  out <- out[do.call(order, ord_args), , drop = FALSE]
  rownames(out) <- NULL
  out
}

aggregate_performance <- function(cell_perf, group_cols, method_order) {
  split_cols <- c(group_cols, "method")
  key <- interaction(cell_perf[split_cols], drop = TRUE, lex.order = TRUE)
  metrics <- c(
    "frob", "frob_sq", "msfe", "rmfe", "hamming", "precision", "recall", "F1",
    "off_hamming", "off_precision", "off_recall", "off_F1",
    "support_hat", "support_hat_offdiag", "rho_hat"
  )
  pieces <- lapply(split(cell_perf, key), function(a) {
    r <- a[1, split_cols, drop = FALSE]
    r$n_cells <- nrow(a)
    for (m in metrics) {
      mean_col <- paste0(m, "_mean"); med_col <- paste0(m, "_median")
      xm <- a[[mean_col]]; xd <- a[[med_col]]
      r[[paste0("avg_cell_mean_", m)]] <- if (any(is.finite(xm))) mean(xm, na.rm = TRUE) else NA_real_
      r[[paste0("avg_cell_median_", m)]] <- if (any(is.finite(xd))) mean(xd, na.rm = TRUE) else NA_real_
    }
    r
  })
  out <- do.call(rbind, pieces)
  ord_args <- c(lapply(group_cols, function(g) out[[g]]), list(match(out$method, method_order)))
  out <- out[do.call(order, ord_args), , drop = FALSE]
  rownames(out) <- NULL
  out
}

sort_aggregate_table <- function(x, cfg, method_order) {
  if (!is.data.frame(x)) stop("sort_aggregate_table(): x must be a data.frame")
  if (!nrow(x)) return(x)
  if (!all(c("tuning", "method") %in% names(x))) stop("tuning and method columns are required")
  ord_args <- list(match(x$tuning, cfg$tunings))
  if ("dgp" %in% names(x)) ord_args[[length(ord_args)+1L]] <- match(x$dgp, cfg$dgps)
  if ("d" %in% names(x)) ord_args[[length(ord_args)+1L]] <- x$d
  if ("T_sample" %in% names(x)) ord_args[[length(ord_args)+1L]] <- x$T_sample
  ord_args[[length(ord_args)+1L]] <- match(x$method, method_order)
  if (any(vapply(ord_args, length, integer(1)) != nrow(x))) stop("Ordering-vector length mismatch")
  x <- x[do.call(order, ord_args), , drop = FALSE]
  rownames(x) <- NULL
  x
}

make_overall_rank_wide <- function(overall, method_order) {
  th <- overall[overall$tuning == "theory", , drop = FALSE]
  cv <- overall[overall$tuning == "cv", , drop = FALSE]
  th <- th[match(method_order, th$method), ]; cv <- cv[match(method_order, cv$method), ]
  data.frame(
    method = method_order,
    VAR_Th_rank_frob = th$avg_rank_frob,
    VAR_Th_rank_rmfe = th$avg_rank_rmfe,
    VAR_Th_rank_hamming = th$avg_rank_hamming,
    VAR_Th_rank_F1 = th$avg_rank_F1,
    VAR_CV_rank_frob = cv$avg_rank_frob,
    VAR_CV_rank_rmfe = cv$avg_rank_rmfe,
    VAR_CV_rank_hamming = cv$avg_rank_hamming,
    VAR_CV_rank_F1 = cv$avg_rank_F1,
    stringsAsFactors = FALSE
  )
}

make_offdiag_rank_wide <- function(overall, method_order) {
  th <- overall[overall$tuning == "theory", , drop = FALSE]
  cv <- overall[overall$tuning == "cv", , drop = FALSE]
  th <- th[match(method_order, th$method), ]; cv <- cv[match(method_order, cv$method), ]
  data.frame(
    method = method_order,
    VAR_Th_rank_off_hamming = th$avg_rank_off_hamming,
    VAR_Th_rank_off_F1 = th$avg_rank_off_F1,
    VAR_CV_rank_off_hamming = cv$avg_rank_off_hamming,
    VAR_CV_rank_off_F1 = cv$avg_rank_off_F1,
    stringsAsFactors = FALSE
  )
}

make_overall_performance_wide <- function(overall, method_order) {
  th <- overall[overall$tuning == "theory", , drop = FALSE]
  cv <- overall[overall$tuning == "cv", , drop = FALSE]
  th <- th[match(method_order, th$method), ]; cv <- cv[match(method_order, cv$method), ]
  core <- c("frob", "rmfe", "hamming", "F1", "off_hamming", "off_F1")
  out <- data.frame(method = method_order, stringsAsFactors = FALSE)
  for (m in core) {
    out[[paste0("VAR_Th_avg_cell_mean_", m)]] <- th[[paste0("avg_cell_mean_", m)]]
    out[[paste0("VAR_Th_avg_cell_median_", m)]] <- th[[paste0("avg_cell_median_", m)]]
    out[[paste0("VAR_CV_avg_cell_mean_", m)]] <- cv[[paste0("avg_cell_mean_", m)]]
    out[[paste0("VAR_CV_avg_cell_median_", m)]] <- cv[[paste0("avg_cell_median_", m)]]
  }
  out
}

write_rank_and_performance_tables <- function(results_dir, cell_perf, cfg, method_order) {
  cell_ranks <- rank_cell_methods(cell_perf, cfg, method_order)
  write_csv(cell_ranks, file.path(results_dir, "ranks", "cell_ranks.csv"))
  write_csv(cell_ranks[cell_ranks$tuning == "theory", , drop = FALSE],
            file.path(results_dir, "ranks", "cell_ranks_theory.csv"))
  write_csv(cell_ranks[cell_ranks$tuning == "cv", , drop = FALSE],
            file.path(results_dir, "ranks", "cell_ranks_cv.csv"))

  rank_overall <- sort_aggregate_table(aggregate_ranks(cell_ranks, c("tuning"), method_order), cfg, method_order)
  rank_dgp <- sort_aggregate_table(aggregate_ranks(cell_ranks, c("tuning", "dgp"), method_order), cfg, method_order)
  rank_d <- sort_aggregate_table(aggregate_ranks(cell_ranks, c("tuning", "d"), method_order), cfg, method_order)
  rank_T <- sort_aggregate_table(aggregate_ranks(cell_ranks, c("tuning", "T_sample"), method_order), cfg, method_order)

  write_csv(rank_overall, file.path(results_dir, "ranks", "average_ranks_overall.csv"))
  write_csv(rank_dgp, file.path(results_dir, "ranks", "average_ranks_by_dgp.csv"))
  write_csv(rank_d, file.path(results_dir, "ranks", "average_ranks_by_d.csv"))
  write_csv(rank_T, file.path(results_dir, "ranks", "average_ranks_by_T.csv"))
  write_csv(make_overall_rank_wide(rank_overall, method_order),
            file.path(results_dir, "ranks", "rank_table_overall_wide.csv"))
  write_csv(make_offdiag_rank_wide(rank_overall, method_order),
            file.path(results_dir, "ranks", "rank_table_offdiag_wide.csv"))

  perf_overall <- sort_aggregate_table(aggregate_performance(cell_perf, c("tuning"), method_order), cfg, method_order)
  perf_dgp <- sort_aggregate_table(aggregate_performance(cell_perf, c("tuning", "dgp"), method_order), cfg, method_order)
  perf_d <- sort_aggregate_table(aggregate_performance(cell_perf, c("tuning", "d"), method_order), cfg, method_order)
  perf_T <- sort_aggregate_table(aggregate_performance(cell_perf, c("tuning", "T_sample"), method_order), cfg, method_order)

  write_csv(perf_overall, file.path(results_dir, "summary", "average_performance_overall.csv"))
  write_csv(perf_dgp, file.path(results_dir, "summary", "average_performance_by_dgp.csv"))
  write_csv(perf_d, file.path(results_dir, "summary", "average_performance_by_d.csv"))
  write_csv(perf_T, file.path(results_dir, "summary", "average_performance_by_T.csv"))
  write_csv(make_overall_performance_wide(perf_overall, method_order),
            file.path(results_dir, "summary", "performance_table_overall_wide.csv"))

  invisible(list(cell_ranks = cell_ranks, rank_overall = rank_overall, perf_overall = perf_overall))
}
