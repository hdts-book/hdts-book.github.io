# Figure 1.1 — Stable versus unit-root bivariate VAR(1)
# High-Dimensional Time Series
#
# DESIGN NOTES
# ------------
# Layout is a single ggplot with facet_grid(component ~ model).  Column strips
# carry the model title, row strips carry the component label.  Strip heights
# are computed by the gtable from the text they contain, so a larger title
# widens the strip band instead of colliding with the panel below.
#
# The y scale is free across rows but shared within a row: at a common vertical
# scale the stable path stays inside a band while the unit-root path wanders.
#
# OUTPUTS
# -------
# figures/fig01_01_var_stability.tex   LaTeX/TikZ source        (USE_TIKZ = TRUE)
# figures/fig01_01_var_stability.pdf   publication-quality figure
# figures/fig01_01_var_stability.png   preview, when conversion is available

# FALSE renders through an ordinary PDF device with plain-text labels, which is
# the fastest way to check that the panels themselves are fine.
USE_TIKZ <- TRUE

# Column-title typography.  base_strip_size is the size used by the row strips
# and by the previous version of the figure; the column title is twice that.
base_strip_size  <- 11
title_size       <- 1.5 * base_strip_size

# The coefficient matrix roughly doubles the length of the title line.  At
# title_size = 20.4 pt the combined line needs about 5.8 in per column, so it is
# off by default and belongs in the LaTeX caption instead.  Set to TRUE only if
# you also raise fig_width below.
SHOW_PHI_IN_TITLE <- FALSE

required_pkgs <- c("ggplot2", "sparseVAR", if (USE_TIKZ) "tikzDevice")
missing_pkgs <- required_pkgs[!vapply(
  required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)
)]
if (length(missing_pkgs) > 0L) {
  stop("Install required package(s): ", paste(missing_pkgs, collapse = ", "))
}

library(ggplot2)

# -----------------------------------------------------------------------------
# 1. Model specification and simulation
# -----------------------------------------------------------------------------
n_obs <- 100L
Sigma_eps <- diag(2)

Phi_stable <- matrix(c(0.7, 0.3,
                       0.1, 0.7), nrow = 2L, byrow = TRUE)

Phi_boundary <- matrix(c(1.0, 0.2,
                         0.0, 0.7), nrow = 2L, byrow = TRUE)

spectral_radius <- function(A) max(Mod(eigen(A, only.values = TRUE)$values))

rho_stable   <- spectral_radius(Phi_stable)
rho_boundary <- spectral_radius(Phi_boundary)

stopifnot(rho_stable < 1, abs(rho_boundary - 1) < 1e-10)

# sparseVAR 0.3.0 contains VAR_sim_cpp(), although it is not yet exported.
# Replace this by sparseVAR::VAR_sim_cpp after it becomes part of the public API.
VAR_sim <- getFromNamespace("VAR_sim_cpp", "sparseVAR")

set.seed(12345)

# A burn-in is meaningful only for the stable process.
X_stable <- VAR_sim(
  T = n_obs,
  A = Phi_stable,
  Sigma = Sigma_eps,
  burn = 200L
)

X_boundary <- VAR_sim(
  T = n_obs,
  A = Phi_boundary,
  Sigma = Sigma_eps,
  burn = 200L
)

# -----------------------------------------------------------------------------
# 2. Strip labels
# -----------------------------------------------------------------------------
# In TikZ mode these strings reach LaTeX verbatim (sanitize = FALSE).
# smallmatrix keeps the coefficient matrix to about one text line.
model_label <- function(name, Phi, rho) {
  if (USE_TIKZ) {
    phi_part <- if (SHOW_PHI_IN_TITLE) {
      sprintf(
        paste0("\\quad $\\Phi_1=\\bigl(\\begin{smallmatrix}",
               "%.1f & %.1f \\\\ %.1f & %.1f",
               "\\end{smallmatrix}\\bigr)$"),
        Phi[1, 1], Phi[1, 2], Phi[2, 1], Phi[2, 2]
      )
    } else ""
    sprintf("%s%s\\quad $\\rho(\\Phi_1)=%.3f$", name, phi_part, rho)
  } else {
    phi_part <- if (SHOW_PHI_IN_TITLE) {
      sprintf("   Phi = (%.1f, %.1f; %.1f, %.1f)",
              Phi[1, 1], Phi[1, 2], Phi[2, 1], Phi[2, 2])
    } else ""
    sprintf("%s%s   rho = %.3f", name, phi_part, rho)
  }
}

lab_stable   <- model_label("Stable VAR(1)",      Phi_stable,   rho_stable)
lab_boundary <- model_label("Unit-root boundary", Phi_boundary, rho_boundary)

comp_labels <- if (USE_TIKZ) c("$X_{1,t}$", "$X_{2,t}$") else c("X1", "X2")
time_label  <- if (USE_TIKZ) "$t$" else "t"

# -----------------------------------------------------------------------------
# 3. Long-format data
# -----------------------------------------------------------------------------
to_long <- function(X, model_lab) {
  data.frame(
    time  = rep(seq_len(ncol(X)), times = nrow(X)),
    value = as.numeric(t(X)),
    comp  = factor(rep(comp_labels, each = ncol(X)), levels = comp_labels),
    model = factor(model_lab, levels = c(lab_stable, lab_boundary))
  )
}

dat <- rbind(
  to_long(X_stable,   lab_stable),
  to_long(X_boundary, lab_boundary)
)

stopifnot(nrow(dat) == 4L * n_obs)

# -----------------------------------------------------------------------------
# 4. The figure
# -----------------------------------------------------------------------------
fig <- ggplot(dat, aes(x = time, y = value)) +
  geom_hline(yintercept = 0, linewidth = 0.28, linetype = "22",
             colour = "grey67") +
  geom_line(linewidth = 0.45, colour = "grey12", lineend = "round") +
  facet_grid(comp ~ model, scales = "free_y") +
  scale_x_continuous(
    breaks = seq(0, n_obs, by = 25),
    expand = expansion(mult = c(0.015, 0.02))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.09, 0.09))) +
  labs(x = time_label, y = NULL) +
  theme_bw(base_size = 10.5, base_family = "serif") +
  theme(
    panel.grid.major = element_line(linewidth = 0.22, colour = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(linewidth = 0.45, colour = "grey30"),
    panel.spacing.x  = unit(9, "pt"),
    panel.spacing.y  = unit(7, "pt"),
    strip.background = element_rect(fill = "grey95", colour = "grey30",
                                    linewidth = 0.45),
    # Column title at twice the row-strip size.  The generous vertical margin
    # keeps the enlarged text off the strip rule.
    strip.text.x     = element_text(size = title_size, face = "bold",
                                    margin = margin(8, 6, 8, 6)),
    strip.text.y     = element_text(size = base_strip_size, angle = 0,
                                    margin = margin(4, 5, 4, 5)),
    axis.text        = element_text(size = 9.2, colour = "grey15"),
    axis.title.x     = element_text(size = 10.6, margin = margin(t = 5)),
    plot.margin      = margin(4, 5, 4, 4)
  )

# -----------------------------------------------------------------------------
# 5. Rendering
# -----------------------------------------------------------------------------
out_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tex_file <- file.path(out_dir, "fig01_01_var_stability.tex")
pdf_file <- file.path(out_dir, "fig01_01_var_stability.pdf")
png_file <- file.path(out_dir, "fig01_01_var_stability.png")

# Taller than before so that the enlarged strip band does not eat into the
# panels: the two strips now take roughly 0.8 in of the total height.
fig_width  <- 9.0
fig_height <- 5.6

if (!USE_TIKZ) {

  grDevices::pdf(pdf_file, width = fig_width, height = fig_height)
  print(fig)
  grDevices::dev.off()
  message("Preview PDF written to ", pdf_file)

} else {

  # smallmatrix and \bigl need amsmath, which is not in the default tikzDevice
  # preamble.  The option is restored afterwards.
  old_pkgs <- getOption("tikzLatexPackages")
  options(tikzLatexPackages = c(
    old_pkgs,
    "\\usepackage{amsmath}\n",
    "\\usepackage{amssymb}\n"
  ))

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

  # An empty figure shows up here as a near-total absence of path commands,
  # long before anyone has to open the PDF.
  tex_src <- readLines(tex_file, warn = FALSE)
  n_paths <- sum(grepl("\\\\path|\\\\draw", tex_src))
  message("TikZ source: ", length(tex_src), " lines, ", n_paths,
          " path commands.")
  if (n_paths < 50L) {
    warning("The TikZ source contains very few drawing commands. ",
            "The panels are probably empty.")
  }

  pdflatex <- Sys.which("pdflatex")
  if (nzchar(pdflatex)) {
    old_wd <- setwd(out_dir)
    log_out <- system2(
      pdflatex,
      args = c("-interaction=nonstopmode", basename(tex_file)),
      stdout = TRUE, stderr = TRUE
    )
    setwd(old_wd)
    if (!file.exists(pdf_file)) {
      stop("LaTeX compilation failed:\n", paste(log_out, collapse = "\n"))
    }
    message("Figure PDF: ", pdf_file)
  } else {
    warning("pdflatex not found; the .tex was written but not compiled.")
  }
}

# Optional PNG preview.
if (file.exists(pdf_file)) {
  if (requireNamespace("pdftools", quietly = TRUE)) {
    pdftools::pdf_convert(
      pdf = pdf_file, format = "png", pages = 1L,
      dpi = 300, filenames = png_file, verbose = FALSE
    )
  } else {
    pdftoppm <- Sys.which("pdftoppm")
    if (nzchar(pdftoppm)) {
      system2(pdftoppm, args = c(
        "-png", "-r", "300", "-singlefile", pdf_file,
        tools::file_path_sans_ext(png_file)
      ))
    } else {
      warning("PNG conversion skipped: install 'pdftools' or 'pdftoppm'.")
    }
  }
}

if (file.exists(png_file)) message("Figure PNG: ", png_file)
