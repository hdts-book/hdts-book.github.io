# General utilities and deterministic seed bookkeeping.

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
  off <- switch(
    dgp,
    dgp1 = 0L,
    dgp2 = 1000L,
    dgp3 = 2000L,
    stop("Unknown DGP: ", dgp)
  )
  as.integer(base + off + idx_d)
}

data_seed <- function(master_seed, dgp, idx_d, idx_T, rep) {
  off <- switch(
    dgp,
    dgp1 = 0L,
    dgp2 = 1000000L,
    dgp3 = 2000000L,
    stop("Unknown DGP: ", dgp)
  )
  as.integer(master_seed + off + 100000L * idx_d + 1000L * idx_T + rep)
}

make_design_manifest <- function(cfg) {
  z <- expand.grid(
    dgp = cfg$dgps,
    idx_d = seq_along(cfg$d_grid),
    idx_T = seq_along(cfg$T_grid),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  z$d <- cfg$d_grid[z$idx_d]
  z$T_sample <- cfg$T_grid[z$idx_T]
  z$truth_seed <- mapply(
    truth_seed, z$dgp, z$idx_d,
    MoreArgs = list(base = cfg$truth_seed_base)
  )
  z <- z[, c("dgp", "idx_d", "d", "idx_T", "T_sample", "truth_seed")]
  z[order(match(z$dgp, cfg$dgps), z$idx_d, z$idx_T), , drop = FALSE]
}

make_seed_manifest <- function(cfg) {
  z <- expand.grid(
    dgp = cfg$dgps,
    idx_d = seq_along(cfg$d_grid),
    idx_T = seq_along(cfg$T_grid),
    rep = seq_len(cfg$nrep),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  z$d <- cfg$d_grid[z$idx_d]
  z$T_sample <- cfg$T_grid[z$idx_T]
  z$truth_seed <- mapply(
    truth_seed, z$dgp, z$idx_d,
    MoreArgs = list(base = cfg$truth_seed_base)
  )
  z$data_seed <- mapply(
    data_seed, z$dgp, z$idx_d, z$idx_T, z$rep,
    MoreArgs = list(master_seed = cfg$master_seed)
  )
  z <- z[, c(
    "dgp", "idx_d", "d", "idx_T", "T_sample", "rep",
    "truth_seed", "data_seed"
  )]
  z <- z[order(match(z$dgp, cfg$dgps), z$idx_d, z$idx_T, z$rep), , drop = FALSE]
  rownames(z) <- NULL
  if (anyDuplicated(z$data_seed)) stop("Duplicate data seeds detected.")
  z
}

require_replication_packages <- function(cfg) {
  required <- c("sparseVAR", "MASS", "igraph", "foreach", "doParallel")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "))

  got <- as.character(utils::packageVersion("sparseVAR"))
  if (isTRUE(cfg$strict_package_version) && got != cfg$required_sparseVAR) {
    stop(
      "This replication requires sparseVAR ", cfg$required_sparseVAR,
      "; installed version is ", got, "."
    )
  }

  ns <- asNamespace("sparseVAR")
  need <- c(
    "sVAR_variants", "sVAR_adalasso_fista", "sVAR_scaled",
    "sVAR_threshold", "sVAR_refit", "sVAR_relax", "sVAR_debias",
    "sVAR_book_control"
  )
  miss <- need[!vapply(need, exists, logical(1), envir = ns, inherits = FALSE)]
  if (length(miss)) stop("sparseVAR 0.4.3 API missing: ", paste(miss, collapse = ", "))

  ctl <- sparseVAR::sVAR_book_control()
  if (!identical(ctl$whiten, "none")) stop("sVAR 0.4.3 book control does not report whiten='none'.")
  if (!identical(ctl$penalize_diag, FALSE)) stop("sVAR 0.4.3 book control does not report penalize_diag=FALSE.")
  invisible(TRUE)
}

safe_reset_results <- function(project_dir, output_dir, clean = TRUE) {
  out <- normalizePath(file.path(project_dir, output_dir), winslash = "/", mustWork = FALSE)
  root <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  if (!startsWith(out, paste0(root, "/"))) stop("Refusing to clean output outside PROJECT_DIR: ", out)
  if (basename(out) != basename(output_dir)) stop("Unexpected output directory: ", out)
  if (dir.exists(out) && isTRUE(clean)) unlink(out, recursive = TRUE, force = TRUE)
  dirs <- c(out, file.path(out, "raw"), file.path(out, "summary"),
            file.path(out, "ranks"), file.path(out, "manifests"))
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  out
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
  out <- lapply(names(cfg), function(nm) {
    v <- cfg[[nm]]
    if (length(v) == 0L) val <- ""
    else if (is.null(v)) val <- "NULL"
    else val <- paste(as.character(v), collapse = "|")
    data.frame(key = nm, value = val, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

package_manifest <- function(cfg) {
  pkgs <- c("sparseVAR", "MASS", "igraph", "foreach", "doParallel", "glmnet", "Matrix")
  versions <- vapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
  }, character(1))
  data.frame(
    item = c("R", "platform", pkgs),
    version = c(as.character(getRversion()), R.version$platform, unname(versions)),
    stringsAsFactors = FALSE
  )
}

code_manifest <- function(project_dir) {
  files <- sort(list.files(project_dir, pattern = "\\.R$", full.names = TRUE))
  data.frame(file = basename(files), md5 = unname(tools::md5sum(files)), stringsAsFactors = FALSE)
}

assert_method_set <- function(methods, method_order) {
  if (!is.list(methods)) stop("sVAR_variants() did not return a method list.")
  if (!setequal(names(methods), method_order)) {
    stop(
      "Unexpected method set. Expected: ", paste(method_order, collapse = ", "),
      "; got: ", paste(names(methods), collapse = ", ")
    )
  }
  invisible(TRUE)
}

assert_production_profile <- function(cfg) {
  if (!identical(cfg$penalize_diag, FALSE)) stop("Production profile must use penalize_diag=FALSE.")
  if (!identical(cfg$whiten, "none")) stop("Production profile must use whiten='none'.")
  if (!isTRUE(cfg$standardize)) stop("Production profile must use standardize=TRUE.")
  if (!identical(cfg$n.cl, 50L)) stop("Production profile must use exactly 50 workers.")
  if (!identical(cfg$threshold_factor, 1.0)) stop("Production profile requires threshold_factor=1.0 (eta=lambda).")
  if (!identical(cfg$rng_version, "4.3.3")) stop("Production profile requires rng_version=\"4.3.3\".")
  invisible(TRUE)
}
