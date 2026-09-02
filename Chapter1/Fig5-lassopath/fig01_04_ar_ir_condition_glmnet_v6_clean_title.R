# Figure 1.4 — LASSO paths and the AR(2)/AR(3) irrepresentability condition
# High-Dimensional Time Series
#
# VERSION 6
# ---------
# - Retains the original plot.glmnet coefficient-path colours and solid lines.
# - Removes the auxiliary IC diagnostic notation from the facet titles.
# - Sets facet titles to 15 pt and both axis titles to 12 pt.
# - Keeps the panels free of IC callouts, arrows, and comment boxes.
# - Uses the original glmnet-style y-axis label, 'Coefficients'.
# - Moves the interpretation formerly shown inside the panels to the LaTeX
#   caption, including the behavior of the inactive third-lag coefficient.
# - Uses a single ggplot object and the same TikZ rendering workflow as the
#   Figure 1.1 reference script.
#
# OUTPUTS
# -------
# figures/fig01_04_ar_ir_condition_glmnet_v6.tex
# figures/fig01_04_ar_ir_condition_glmnet_v6.pdf
# figures/fig01_04_ar_ir_condition_glmnet_v6.png
# figures/fig01_04_ar_ir_condition_caption_v6.tex
# figures/fig01_04_ar_ir_condition_book_insert_v6.tex

USE_TIKZ <- TRUE

required_pkgs <- c("ggplot2", "glmnet", if (USE_TIKZ) "tikzDevice")
missing_pkgs <- required_pkgs[!vapply(
  required_pkgs,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_pkgs) > 0L) {
  stop(
    "Install required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    call. = FALSE
  )
}

library(ggplot2)

# -----------------------------------------------------------------------------
# 1. Simulation and the original AR(3) design
# -----------------------------------------------------------------------------
n_obs           <- 500L
p_fit           <- 3L
innovation_sd   <- 1
simulation_seed <- 1234L

cases <- list(
  list(
    case_id = "satisfied",
    status  = "IC satisfied",
    phi     = c(0.2, 0.5)
  ),
  list(
    case_id = "violated",
    status  = "IC violated",
    phi     = c(-1.0, -0.6)
  )
)

# AR(2) causality region used in the chapter:
# |phi_2| < 1, phi_2 + phi_1 < 1, phi_2 - phi_1 < 1.
is_causal_ar2 <- function(phi) {
  stopifnot(length(phi) == 2L)
  abs(phi[2L]) < 1 &&
    phi[2L] + phi[1L] < 1 &&
    phi[2L] - phi[1L] < 1
}

# For an AR(2) truth embedded in an AR(3) fit, the chapter's IC quantity is
# |phi_2 sign(phi_1) + phi_1 sign(phi_2)|.
ic_value_ar2_in_ar3 <- function(phi) {
  stopifnot(length(phi) == 2L)
  abs(phi[2L] * sign(phi[1L]) + phi[1L] * sign(phi[2L]))
}

# Reproduce the lag construction in the original code: unavailable presample
# lags are filled with the sample mean, so X and y both contain n observations.
build_original_lag_design <- function(y, p) {
  y  <- as.numeric(y)
  n  <- length(y)
  mu <- mean(y)
  id <- seq_len(n)

  X <- vapply(seq_len(p), function(j) {
    lag_id <- id - j
    c(rep(mu, sum(lag_id <= 0L)), y[lag_id[lag_id > 0L]])
  }, FUN.VALUE = numeric(n))

  colnames(X) <- paste0("lag", seq_len(p))
  list(x = X, y = y)
}

lag_ids <- paste0("lag", seq_len(p_fit))

if (USE_TIKZ) {
  lag_labels <- c(
    "$\\widehat{\\phi}_1$",
    "$\\widehat{\\phi}_2$",
    "$\\widehat{\\phi}_3$"
  )
  x_axis_label <- "$\\lVert\\widehat{\\boldsymbol{\\phi}}(\\lambda)\\rVert_1$"
  y_axis_label <- "Coefficients"
} else {
  lag_labels   <- c("phi1", "phi2", "phi3")
  x_axis_label <- "L1 norm"
  y_axis_label <- "Coefficients"
}
names(lag_labels) <- lag_ids

make_case_label <- function(status, phi) {
  if (USE_TIKZ) {
    sprintf(
      "%s\\quad $(\\phi_1,\\phi_2)=(%.1f,%.1f)$",
      status, phi[1L], phi[2L]
    )
  } else {
    sprintf(
      "%s   (phi1, phi2) = (%.1f, %.1f)",
      status, phi[1L], phi[2L]
    )
  }
}

fit_one_case <- function(case_spec) {
  phi      <- case_spec$phi
  ic_value <- ic_value_ar2_in_ar3(phi)

  if (!is_causal_ar2(phi)) {
    stop(
      "The supplied AR(2) parameter is not causal: ",
      paste(phi, collapse = ", "),
      call. = FALSE
    )
  }

  set.seed(simulation_seed)
  y <- stats::arima.sim(
    n = n_obs,
    model = list(ar = phi),
    sd = innovation_sd
  )

  design <- build_original_lag_design(y, p = p_fit)

  fit <- glmnet::glmnet(
    x = design$x,
    y = design$y,
    family = "gaussian",
    alpha = 1,
    intercept = TRUE,
    standardize = TRUE,
    nlambda = 100L
  )

  beta <- as.matrix(fit$beta)
  if (nrow(beta) != p_fit) {
    stop("Unexpected number of glmnet coefficient paths.", call. = FALSE)
  }

  # Match plot.glmnet(xvar = "norm"): horizontal coordinate is the L1 norm
  # of the fitted coefficient vector at each value of lambda.
  l1_norm <- colSums(abs(beta))
  n_path  <- length(l1_norm)

  case_label <- make_case_label(case_spec$status, phi)

  path_data <- data.frame(
    case_id    = case_spec$case_id,
    case_label = case_label,
    path_index = rep(seq_len(n_path), times = p_fit),
    l1_norm    = rep(l1_norm, times = p_fit),
    estimate   = as.vector(t(beta)),
    lag_id     = factor(rep(lag_ids, each = n_path), levels = lag_ids),
    stringsAsFactors = FALSE
  )

  diagnostics <- data.frame(
    case_id   = case_spec$case_id,
    phi1      = phi[1L],
    phi2      = phi[2L],
    causal_CC = TRUE,
    IC_value  = ic_value,
    IC_holds  = ic_value < 1,
    stringsAsFactors = FALSE
  )

  list(path = path_data, diagnostics = diagnostics, fit = fit)
}

results     <- lapply(cases, fit_one_case)
path_dat    <- do.call(rbind, lapply(results, `[[`, "path"))
diagnostics <- do.call(rbind, lapply(results, `[[`, "diagnostics"))

case_levels <- vapply(cases, function(z) {
  make_case_label(z$status, z$phi)
}, FUN.VALUE = character(1))
path_dat$case_label <- factor(path_dat$case_label, levels = case_levels)

stopifnot(
  diagnostics$IC_holds[diagnostics$case_id == "satisfied"],
  !diagnostics$IC_holds[diagnostics$case_id == "violated"],
  all(is.finite(path_dat$l1_norm)),
  all(is.finite(path_dat$estimate))
)

print(diagnostics, row.names = FALSE)

# -----------------------------------------------------------------------------
# 2. Original glmnet colours and direct path labels
# -----------------------------------------------------------------------------
# plot.glmnet() calls matplot(..., lty = 1) without supplying 'col'.  matplot
# therefore uses integer colours 1, 2, 3, i.e. the first three entries of the
# active R palette. Reading palette() reproduces the colours in the same R
# session.
glmnet_colours <- grDevices::palette()[seq_len(p_fit)]
names(glmnet_colours) <- lag_ids

# Label each coefficient at the right-hand endpoint of its path.
endpoint_rows <- unlist(lapply(
  split(
    seq_len(nrow(path_dat)),
    interaction(path_dat$case_id, path_dat$lag_id, drop = TRUE)
  ),
  function(ii) ii[which.max(path_dat$l1_norm[ii])]
), use.names = FALSE)

endpoint_dat <- path_dat[endpoint_rows, , drop = FALSE]
endpoint_dat$path_label <- unname(lag_labels[as.character(endpoint_dat$lag_id)])

# -----------------------------------------------------------------------------
# 3. ggplot2 coefficient-path figure
# -----------------------------------------------------------------------------
fig <- ggplot(
  path_dat,
  aes(
    x = l1_norm,
    y = estimate,
    group = lag_id,
    colour = lag_id
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.30,
    linetype = "22",
    colour = "grey70"
  ) +
  # All paths are solid, matching plot.glmnet(), which uses lty = 1.
  geom_path(linewidth = 0.88, lineend = "round") +
  geom_text(
    data = endpoint_dat,
    aes(
      x = l1_norm,
      y = estimate,
      label = path_label,
      colour = lag_id
    ),
    inherit.aes = FALSE,
    hjust = -0.12,
    vjust = -0.25,
    size = 3.55,
    family = "serif",
    show.legend = FALSE
  ) +
  facet_wrap(
    vars(case_label),
    nrow = 1L,
    scales = "free"
  ) +
  scale_colour_manual(
    values = glmnet_colours,
    breaks = lag_ids,
    labels = unname(lag_labels),
    drop = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.018, 0.18))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.10, 0.16))
  ) +
  labs(x = x_axis_label, y = y_axis_label) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 10.5, base_family = "serif") +
  theme(
    panel.grid.major = element_line(linewidth = 0.22, colour = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.48, colour = "grey30"),
    panel.spacing.x = grid::unit(12, "pt"),
    strip.background = element_rect(
      fill = "grey95",
      colour = "grey30",
      linewidth = 0.48
    ),
    strip.text.x = element_text(
      size = 15.0,
      face = "bold",
      margin = margin(8, 7, 8, 7)
    ),
    axis.text = element_text(size = 10.0, colour = "grey15"),
    axis.title.x = element_text(
      size = 12.0,
      margin = margin(t = 6)
    ),
    # Conventional glmnet-style y-axis title.
    axis.title.y = element_text(
      size = 12.0,
      angle = 90,
      vjust = 0.5,
      margin = margin(r = 7)
    ),
    legend.position = "none",
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank(),
    plot.margin = margin(5, 16, 5, 8)
  )

# -----------------------------------------------------------------------------
# 4. Rendering
# -----------------------------------------------------------------------------
out_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

base_name <- "fig01_04_ar_ir_condition_glmnet_v6"
tex_file  <- file.path(out_dir, paste0(base_name, ".tex"))
pdf_file  <- file.path(out_dir, paste0(base_name, ".pdf"))
png_file  <- file.path(out_dir, paste0(base_name, ".png"))
caption_file <- file.path(
  out_dir,
  "fig01_04_ar_ir_condition_caption_v6.tex"
)
insert_file <- file.path(
  out_dir,
  "fig01_04_ar_ir_condition_book_insert_v6.tex"
)

# Remove stale outputs before rendering.
unlink(c(
  tex_file,
  pdf_file,
  png_file,
  caption_file,
  insert_file,
  file.path(out_dir, paste0(base_name, ".aux")),
  file.path(out_dir, paste0(base_name, ".log"))
), force = TRUE)

# Keep the interpretation out of the plotting panels and place it in the book
# caption. The first file contains only \caption and \label; the second is a
# complete figure environment ready to insert into the chapter source.
caption_lines <- c(
  "\\caption{LASSO coefficient paths for an AR$(2)$ process fitted by an AR$(3)$ model. ",
  "Left: $(\\phi_1,\\phi_2)=(0.2,0.5)$, which satisfies both the causality condition (CC) ",
  "and the irrepresentability condition (IC). Right: $(\\phi_1,\\phi_2)=(-1,-0.6)$, ",
  "which is causal but violates IC. The horizontal axis is the $\\ell_1$ norm along ",
  "the solution path, and the three curves correspond to the three fitted lag ",
  "coefficients. In both fits the true third-lag coefficient is zero. When IC holds ",
  "(left), $\\widehat{\\phi}_3$ remains close to zero along the path; when IC fails ",
  "(right), the inactive third-lag coefficient enters the LASSO path, illustrating ",
  "that causality alone does not guarantee exclusion of a spurious lag.}",
  "\\label{fig:vars-ar2-lasso-paths}"
)
writeLines(caption_lines, caption_file)

insert_lines <- c(
  "\\begin{figure}[tbp]",
  "  \\centering",
  sprintf(
    "  \\includegraphics[width=\\textwidth]{figures/%s.pdf}",
    base_name
  ),
  paste0("  ", caption_lines),
  "\\end{figure}"
)
writeLines(insert_lines, insert_file)

# Slightly wider than v2 to accommodate the horizontal y-axis title cleanly.
fig_width  <- 9.3
fig_height <- 4.55

compile_standalone_tex <- function(tex_path, pdf_path) {
  pdflatex <- Sys.which("pdflatex")
  if (!nzchar(pdflatex)) {
    stop(
      "pdflatex was not found. Install a LaTeX distribution and rerun the script.",
      call. = FALSE
    )
  }

  original_wd <- getwd()
  on.exit(setwd(original_wd), add = TRUE)
  setwd(dirname(tex_path))

  latex_output <- system2(
    pdflatex,
    args = c(
      "-interaction=nonstopmode",
      "-halt-on-error",
      basename(tex_path)
    ),
    stdout = TRUE,
    stderr = TRUE
  )

  if (!file.exists(basename(pdf_path))) {
    stop(
      "LaTeX compilation failed:\n",
      paste(latex_output, collapse = "\n"),
      call. = FALSE
    )
  }

  invisible(latex_output)
}

if (!USE_TIKZ) {
  grDevices::pdf(pdf_file, width = fig_width, height = fig_height)
  print(fig)
  grDevices::dev.off()
} else {
  old_pkgs <- getOption("tikzLatexPackages")
  options(tikzLatexPackages = unique(c(
    old_pkgs,
    "\\usepackage{amsmath}\n",
    "\\usepackage{amssymb}\n"
  )))

  tikzDevice::tikz(
    file = tex_file,
    width = fig_width,
    height = fig_height,
    standAlone = TRUE,
    sanitize = FALSE,
    engine = "pdftex"
  )
  print(fig)
  grDevices::dev.off()

  options(tikzLatexPackages = old_pkgs)

  tex_src <- readLines(tex_file, warn = FALSE)
  n_paths <- sum(grepl("\\\\path|\\\\draw", tex_src))
  message("TikZ source: ", length(tex_src), " lines, ", n_paths, " path commands.")
  if (n_paths < 40L) {
    warning("The TikZ source contains unexpectedly few drawing commands.")
  }

  compile_standalone_tex(tex_file, pdf_file)
}

# Optional 300-dpi PNG for GitHub Pages.
if (file.exists(pdf_file)) {
  if (requireNamespace("pdftools", quietly = TRUE)) {
    pdftools::pdf_convert(
      pdf = pdf_file,
      format = "png",
      pages = 1L,
      dpi = 300,
      filenames = png_file,
      verbose = FALSE
    )
  } else {
    pdftoppm <- Sys.which("pdftoppm")
    if (nzchar(pdftoppm)) {
      system2(
        pdftoppm,
        args = c(
          "-png", "-r", "300", "-singlefile",
          pdf_file,
          tools::file_path_sans_ext(png_file)
        )
      )
    } else {
      warning("PNG conversion skipped: install 'pdftools' or 'pdftoppm'.")
    }
  }
}

message("DONE — Figure 1.4 v6 generated.")
if (USE_TIKZ) message("TikZ: ", tex_file)
message("PDF:  ", pdf_file)
if (file.exists(png_file)) message("PNG:  ", png_file)
message("Caption: ", caption_file)
message("Book insertion: ", insert_file)
