# sparseVAR 0.4.3 installation + smoke test (portable)

pkg_path <- file.path(getwd(), "sparseVAR_0.4.3_fixed.tar.gz")
stopifnot(file.exists(pkg_path))
install.packages(pkg_path, repos = NULL, type = "source")

library(sparseVAR)
cat("Installed sparseVAR version:", as.character(packageVersion("sparseVAR")), "\n")
stopifnot(as.character(packageVersion("sparseVAR")) == "0.4.3")

ctrl <- sVAR_book_control()
print(ctrl)
stopifnot(identical(ctrl$whiten, "none"),
          identical(ctrl$penalize_diag, FALSE),
          identical(ctrl$standardize, TRUE))

RNGversion("4.3.3")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20260917)
Y <- matrix(rnorm(4 * 120), nrow = 4)
fit <- sVAR_variants(Yt = Y, p = 1, tuning = "theory",
                     penalize_diag = FALSE, standardize = TRUE, whiten = "none",
                     threshold_factor = 1, nlambda = 20L, max_iter = 300L, tol = 1e-6)
stopifnot(identical(fit$control$whiten, "none"),
          identical(fit$control$penalize_diag, FALSE),
          length(fit$methods) == 11L,
          all(vapply(fit$methods, function(z) is.matrix(z) && all(is.finite(z)), logical(1))))
cat("\nSmoke test: PASS\n")
