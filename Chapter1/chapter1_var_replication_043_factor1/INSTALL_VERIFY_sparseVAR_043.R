# sparseVAR 0.4.3 installation + smoke test for the public replication bundle.

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit)) return(dirname(normalizePath(sub("^--file=", "", hit[1]), winslash = "/")))
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) return(dirname(normalizePath(ofile, winslash = "/")))
  normalizePath(getwd(), winslash = "/")
}

PROJECT_DIR <- script_dir()
pkg_path <- file.path(PROJECT_DIR, "sparseVAR_0.4.3_fixed.tar.gz")
stopifnot(file.exists(pkg_path))

install.packages(pkg_path, repos = NULL, type = "source")
suppressPackageStartupMessages(library(sparseVAR))

cat("Installed sparseVAR version: ", as.character(packageVersion("sparseVAR")), "\n", sep = "")
stopifnot(as.character(packageVersion("sparseVAR")) == "0.4.3")

ctrl <- sVAR_book_control()
print(ctrl)
stopifnot(
  identical(ctrl$whiten, "none"),
  identical(ctrl$penalize_diag, FALSE),
  identical(ctrl$standardize, TRUE),
  identical(ctrl$threshold_factor, 1)
)

RNGversion("4.3.3")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20260917)
Y <- matrix(rnorm(4 * 120), nrow = 4)

fit <- sVAR_variants(
  Yt = Y, p = 1, tuning = "theory",
  penalize_diag = FALSE, standardize = TRUE, whiten = "none",
  threshold_factor = 1,
  nlambda = 20L, max_iter = 300L, tol = 1e-6
)

stopifnot(
  identical(fit$control$whiten, "none"),
  identical(fit$control$penalize_diag, FALSE),
  identical(fit$control$threshold_factor, 1),
  length(fit$methods) == 11L,
  all(vapply(fit$methods, function(z) is.matrix(z) && all(is.finite(z)), logical(1)))
)

cat("\nSmoke test: PASS\n")
cat("Methods: ", paste(names(fit$methods), collapse = ", "), "\n", sep = "")
