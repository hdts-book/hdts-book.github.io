# Average ranks within Lasso and Adaptive families, plus winner-map data.

rank_algorithm_cells <- function(cell, cfg, solver_order, family_order) {
  key <- interaction(cell[c("tuning","family","dgp","d","T_sample")], drop=TRUE, lex.order=TRUE)
  out <- do.call(rbind, lapply(split(cell,key), function(a) {
    if (!setequal(a$solver, solver_order)) stop("Incomplete solver set.")
    a <- a[match(solver_order,a$solver),,drop=FALSE]
    # The matched-objective comparison should produce the same estimator up to
    # numerical tolerance.  Round only for ranking/winner assignment so machine-
    # precision noise is not misreported as a statistical solver advantage.
    a$rank_frob <- rank(signif(a$frob_median, 8), ties.method="average")
    a$rank_rmfe <- rank(signif(a$rmfe_median, 8), ties.method="average")
    a$rank_hamming <- rank(signif(a$hamming_median, 8), ties.method="average")
    a$rank_F1 <- rank(-signif(a$F1_median, 8), ties.method="average")
    a
  }))
  rownames(out) <- NULL
  out
}

aggregate_algorithm_ranks <- function(rr, cfg, solver_order, family_order) {
  key <- interaction(rr[c("tuning","family","solver")], drop=TRUE, lex.order=TRUE)
  out <- do.call(rbind, lapply(split(rr,key), function(a) {
    data.frame(
      tuning=a$tuning[1], family=a$family[1], solver=a$solver[1], n_cells=nrow(a),
      avg_rank_frob=mean(a$rank_frob), avg_rank_rmfe=mean(a$rank_rmfe),
      avg_rank_hamming=mean(a$rank_hamming), avg_rank_F1=mean(a$rank_F1),
      stringsAsFactors=FALSE)
  }))
  out$aggregate <- rowMeans(out[,c("avg_rank_frob","avg_rank_rmfe","avg_rank_hamming","avg_rank_F1")])
  out <- out[order(match(out$tuning,cfg$tunings), match(out$family,family_order), match(out$solver,solver_order)),]
  rownames(out) <- NULL
  out
}

make_winner_map_data <- function(cell, cfg, solver_order, family_order) {
  metrics <- c("frob","rmfe","hamming","F1")
  pieces <- list(); kk <- 1L
  key <- interaction(cell[c("tuning","family","dgp","d","T_sample")], drop=TRUE, lex.order=TRUE)
  for (a in split(cell,key)) {
    for (m in metrics) {
      x <- a[[paste0(m,"_median")]]
      target <- if (m == "F1") max(x) else min(x)
      tol <- 1e-8 * max(1, abs(target))
      win <- a$solver[abs(x-target) <= tol]
      pieces[[kk]] <- data.frame(
        tuning=a$tuning[1], family=a$family[1], dgp=a$dgp[1], d=a$d[1], T_sample=a$T_sample[1],
        criterion=m, winner=paste(win,collapse="+"), n_tied=length(win),
        stringsAsFactors=FALSE)
      kk <- kk+1L
    }
  }
  do.call(rbind,pieces)
}

write_algorithm_tables <- function(results_dir, cell, cfg, solver_order, family_order) {
  rr <- rank_algorithm_cells(cell,cfg,solver_order,family_order)
  avg <- aggregate_algorithm_ranks(rr,cfg,solver_order,family_order)
  winners <- make_winner_map_data(cell,cfg,solver_order,family_order)
  write_csv(rr,file.path(results_dir,"ranks","cell_ranks.csv"))
  write_csv(avg,file.path(results_dir,"ranks","average_ranks.csv"))
  write_csv(avg[avg$tuning=="theory",],file.path(results_dir,"ranks","rank_table_theory.csv"))
  write_csv(avg[avg$tuning=="cv",],file.path(results_dir,"ranks","rank_table_cv.csv"))
  write_csv(winners,file.path(results_dir,"summary","winner_map_data.csv"))

  # Numerical agreement summary: lower is better for all columns.
  key <- interaction(cell[c("tuning","family","solver")],drop=TRUE,lex.order=TRUE)
  num <- do.call(rbind,lapply(split(cell,key),function(a) data.frame(
    tuning=a$tuning[1],family=a$family[1],solver=a$solver[1],
    avg_objective_gap=mean(a$objective_gap_mean,na.rm=TRUE),
    avg_kkt=mean(a$kkt_max_mean,na.rm=TRUE),
    avg_B_diff_to_fista=mean(a$B_diff_to_fista_mean,na.rm=TRUE),
    avg_runtime_seconds=mean(a$runtime_seconds_mean,na.rm=TRUE),
    admm_convergence_rate=if(a$solver[1]=="admm") mean(a$admm_convergence_rate,na.rm=TRUE) else NA_real_,
    stringsAsFactors=FALSE)))
  write_csv(num,file.path(results_dir,"summary","numerical_agreement.csv"))
  invisible(list(cell_ranks=rr,average_ranks=avg,winners=winners,numerical=num))
}
