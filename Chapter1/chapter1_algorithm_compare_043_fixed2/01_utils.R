# Reproducibility, manifests, and common helpers.

`%||%` <- function(x, y) if (is.null(x)) y else x

set_rng_seed <- function(seed, cfg) {
  RNGversion(cfg$rng_version)
  RNGkind(
    kind = cfg$rng_kind,
    normal.kind = cfg$rng_normal_kind,
    sample.kind = cfg$rng_sample_kind
  )
  set.seed(as.integer(seed))
  invisible(seed)
}

largest_singular_value <- function(A) max(svd(A, nu = 0, nv = 0)$d)

spectral_radius <- function(A) {
  if (any(!is.finite(A))) return(NA_real_)
  max(Mod(eigen(A, only.values = TRUE)$values))
}

truth_seed <- function(dgp, idx_d, base) {
  off <- switch(dgp, dgp1 = 0L, dgp2 = 1000L, dgp3 = 2000L,
                stop("Unknown DGP: ", dgp))
  as.integer(base + off + idx_d)
}

data_seed <- function(master_seed, dgp, idx_d, idx_T, rep) {
  off <- switch(dgp, dgp1 = 0L, dgp2 = 1000000L, dgp3 = 2000000L,
                stop("Unknown DGP: ", dgp))
  as.integer(master_seed + off + 100000L * idx_d + 1000L * idx_T + rep)
}

make_design_manifest <- function(cfg) {
  z <- expand.grid(
    dgp = cfg$dgps,
    idx_d = seq_along(cfg$d_grid),
    idx_T = seq_along(cfg$T_grid),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  z$d <- cfg$d_grid[z$idx_d]
  z$T_sample <- cfg$T_grid[z$idx_T]
  z$truth_seed <- mapply(truth_seed, z$dgp, z$idx_d,
                         MoreArgs = list(base = cfg$truth_seed_base))
  z <- z[order(match(z$dgp, cfg$dgps), z$idx_d, z$idx_T), , drop = FALSE]
  rownames(z) <- NULL
  z
}

make_seed_manifest <- function(cfg) {
  z <- expand.grid(
    dgp = cfg$dgps,
    idx_d = seq_along(cfg$d_grid),
    idx_T = seq_along(cfg$T_grid),
    rep = seq_len(cfg$nrep),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  z$d <- cfg$d_grid[z$idx_d]
  z$T_sample <- cfg$T_grid[z$idx_T]
  z$truth_seed <- mapply(truth_seed, z$dgp, z$idx_d,
                         MoreArgs = list(base = cfg$truth_seed_base))
  z$data_seed <- mapply(data_seed, z$dgp, z$idx_d, z$idx_T, z$rep,
                        MoreArgs = list(master_seed = cfg$master_seed))
  z <- z[order(match(z$dgp, cfg$dgps), z$idx_d, z$idx_T, z$rep), , drop = FALSE]
  rownames(z) <- NULL
  if (anyDuplicated(z$data_seed)) stop("Duplicate data seeds detected.")
  z
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "NA")
  invisible(path)
}

append_csv <- function(x, path, write_header = FALSE) {
  utils::write.table(
    x, file = path, sep = ",", row.names = FALSE,
    col.names = write_header, append = !write_header,
    quote = TRUE, na = "NA", qmethod = "double"
  )
  invisible(path)
}

flatten_config <- function(cfg) {
  do.call(rbind, lapply(names(cfg), function(nm) {
    v <- cfg[[nm]]
    val <- if (is.null(v)) "NULL" else paste(as.character(v), collapse = "|")
    data.frame(key = nm, value = val, stringsAsFactors = FALSE)
  }))
}

package_manifest <- function() {
  pkgs <- c("sparseVAR", "glmnet", "MASS", "igraph", "foreach", "doParallel", "Matrix")
  vv <- vapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
  }, character(1))
  data.frame(
    item = c("R", "platform", pkgs),
    version = c(as.character(getRversion()), R.version$platform, unname(vv)),
    stringsAsFactors = FALSE
  )
}

code_manifest <- function(project_dir) {
  ff <- sort(list.files(project_dir, pattern = "\\.[Rr]$", full.names = TRUE))
  data.frame(file = basename(ff), md5 = unname(tools::md5sum(ff)), stringsAsFactors = FALSE)
}

safe_reset_results <- function(project_dir, output_dir, clean = TRUE) {
  out <- normalizePath(file.path(project_dir, output_dir), winslash = "/", mustWork = FALSE)
  root <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  if (!startsWith(out, paste0(root, "/"))) stop("Refusing to clean output outside PROJECT_DIR: ", out)
  if (dir.exists(out) && isTRUE(clean)) unlink(out, recursive = TRUE, force = TRUE)
  for (dd in c("", "raw", "summary", "ranks", "manifests")) {
    dir.create(file.path(out, dd), recursive = TRUE, showWarnings = FALSE)
  }
  out
}

require_algorithm_packages <- function(cfg) {
  req <- c("sparseVAR", "glmnet", "MASS", "igraph", "foreach", "doParallel", "Matrix")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing package(s): ", paste(miss, collapse = ", "))
  got <- as.character(utils::packageVersion("sparseVAR"))
  if (isTRUE(cfg$strict_package_version) && got != cfg$required_sparseVAR) {
    stop("Required sparseVAR ", cfg$required_sparseVAR, "; installed ", got)
  }
  ns <- asNamespace("sparseVAR")
  need <- c("sVAR_adalasso_fista", "fista_lasso_multi_cpp", "admm_varp_multi_cpp", "var_design")
  miss2 <- need[!vapply(need, exists, logical(1), envir = ns, inherits = FALSE)]
  if (length(miss2)) stop("sparseVAR 0.4.3 API/backend missing: ", paste(miss2, collapse = ", "))
  invisible(TRUE)
}
