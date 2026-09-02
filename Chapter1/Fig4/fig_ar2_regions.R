# AR(2) coefficient regions — CC and CC + IC
#
# The stationarity (causality) triangle has vertices (-2,-1), (2,-1), (0,1).
# The inner region is the l1 ball |phi_1| + |phi_2| < 1, a diamond whose upper
# two edges coincide with the upper edges of the triangle.
#
# SIZING
# ------
# tikzDevice typesets the figure at the device size, so the device must already
# be the final printed size.  Rendering large and then shrinking with
# \includegraphics scales the type down with it: a 9 pt label drawn on an 8 in
# canvas and scaled to 3 in prints at about 3.4 pt.
#
# Set TEXTWIDTH_IN to the book's \the\textwidth (put \the\textwidth in the
# source and read it off the log), then include the figure at its natural size:
#
#   \includegraphics[width=0.75\textwidth]{figures/fig_ar2_regions.pdf}
#
# which is a no-op rescale because the widths agree.

USE_TIKZ <- TRUE

TEXTWIDTH_IN <- 4.80      # book \textwidth in inches; 4.80 in ~= 347 pt
WIDTH_FRAC   <- 0.75

required_pkgs <- c("ggplot2", if (USE_TIKZ) "tikzDevice")
missing_pkgs <- required_pkgs[!vapply(
  required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)
)]
if (length(missing_pkgs) > 0L) {
  stop("Install required package(s): ", paste(missing_pkgs, collapse = ", "))
}

library(ggplot2)

# -----------------------------------------------------------------------------
# 1. Typography for the target width
# -----------------------------------------------------------------------------
# At 0.75 \textwidth the drawing area is roughly 3.6 in.  These sizes were
# chosen against that width; scale them together if WIDTH_FRAC changes.
size_title  <- 11.0
size_axis   <-  9.0
size_tick   <-  8.0
size_legend <-  8.5
size_region <-  2.85   # geom_text is in mm: 2.85 mm ~= 8 pt

# -----------------------------------------------------------------------------
# 2. Labels
# -----------------------------------------------------------------------------
# The condition is carried in the legend itself, which is why the figure is set
# at 0.75 \textwidth: at 0.65 the legend row fits but sits uncomfortably close
# to the x axis title.
if (USE_TIKZ) {
  lab_inner <- "CC $+$ IC: $|\\phi_1| + |\\phi_2| < 1$"
  lab_outer <- "CC only"
  lab_x     <- "$\\phi_1$"
  lab_y     <- "$\\phi_2$"
  lab_title <- "Stationarity regions, $p = 2$"
  txt_inner <- "CC $+$ IC"
  txt_outer <- "CC"
} else {
  lab_inner <- "CC + IC:  |phi_1| + |phi_2| < 1"
  lab_outer <- "CC only"
  lab_x     <- "phi_1"
  lab_y     <- "phi_2"
  lab_title <- "Stationarity regions, p = 2"
  txt_inner <- "CC + IC"
  txt_outer <- "CC"
}

region_levels <- c(lab_inner, lab_outer)

# -----------------------------------------------------------------------------
# 3. Palette
# -----------------------------------------------------------------------------
# A: sequential blues.  The nesting of the two regions is carried by lightness,
#    so the distinction survives greyscale printing and colour-vision deficiency.
# B: sage and apricot.  Brighter, but the two fills are close in lightness and
#    blur together in a black-and-white print.
pal_blue <- c(inner = "#6BAED6", outer = "#DEEBF7", line = "#1F3B57")
pal_warm <- c(inner = "#7FC4B0", outer = "#FDE3C0", line = "#34504B")

pal <- pal_blue

# -----------------------------------------------------------------------------
# 4. Region polygons
# -----------------------------------------------------------------------------
poly <- function(x, y, id, region) {
  data.frame(phi1 = x, phi2 = y, id = id,
             region = factor(region, levels = region_levels))
}

regions <- rbind(
  # l1 ball
  poly(c( 0,  1,  0, -1), c( 1,  0, -1,  0), "diamond", lab_inner),
  # stationarity triangle minus the l1 ball
  poly(c(-2, -1,  0),     c(-1,  0, -1),     "left",    lab_outer),
  poly(c( 2,  1,  0),     c(-1,  0, -1),     "right",   lab_outer)
)

# Boundaries drawn as paths so that line type carries meaning: solid for the
# stationarity boundary, dashed for the inner one.
bound_outer <- data.frame(phi1 = c(-2, 0, 2, -2),
                          phi2 = c(-1, 1, -1, -1))

bound_inner <- data.frame(phi1 = c(-1,  0, 1),
                          phi2 = c( 0, -1, 0))

region_text <- data.frame(
  phi1  = c(0, -1.15, 1.15),
  phi2  = c(-0.20, -0.66, -0.66),
  label = c(txt_inner, txt_outer, txt_outer)
)

# -----------------------------------------------------------------------------
# 5. Figure
# -----------------------------------------------------------------------------
fig <- ggplot() +
  geom_polygon(
    data = regions,
    aes(x = phi1, y = phi2, group = id, fill = region),
    colour = NA
  ) +
  geom_path(
    data = bound_outer, aes(x = phi1, y = phi2),
    linewidth = 0.55, colour = pal[["line"]]
  ) +
  geom_path(
    data = bound_inner, aes(x = phi1, y = phi2),
    linewidth = 0.55, colour = pal[["line"]], linetype = "22"
  ) +
  geom_text(
    data = region_text, aes(x = phi1, y = phi2, label = label),
    size = size_region, colour = pal[["line"]], family = "serif"
  ) +
  scale_fill_manual(
    values = setNames(c(pal[["inner"]], pal[["outer"]]), region_levels),
    breaks = region_levels,
    name = NULL
  ) +
  # Integer breaks on the wide axis, half steps on the short one.
  scale_x_continuous(breaks = seq(-2, 2, by = 1), expand = expansion(0)) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5), expand = expansion(0)) +
  # Honest geometry: MATLAB's default aspect ratio distorts the triangle.
  coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-1, 1), clip = "off") +
  labs(x = lab_x, y = lab_y, title = lab_title) +
  theme_bw(base_size = size_axis, base_family = "serif") +
  theme(
    panel.grid.major = element_line(linewidth = 0.18, colour = "grey87"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(linewidth = 0.45, colour = "grey45"),
    plot.title       = element_text(size = size_title, face = "bold",
                                    hjust = 0.5, margin = margin(b = 5)),
    axis.title.x     = element_text(size = size_axis, margin = margin(t = 3)),
    # Upright reading direction for the vertical axis.
    axis.title.y     = element_text(size = size_axis, angle = 0,
                                    vjust = .5, margin = margin(r = 3)),
    axis.text        = element_text(size = size_tick, colour = "grey20"),
    legend.position  = "bottom",
    legend.text      = element_text(size = size_legend),
    legend.key.size  = unit(8, "pt"),
    legend.spacing.x = unit(4, "pt"),
    legend.margin    = margin(t = 0, b = 0),
    legend.box.spacing = unit(5, "pt"),
    plot.margin      = margin(3, 4, 2, 3)
  )

# -----------------------------------------------------------------------------
# 6. Rendering
# -----------------------------------------------------------------------------
out_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tex_file <- file.path(out_dir, "fig_ar2_regions.tex")
pdf_file <- file.path(out_dir, "fig_ar2_regions.pdf")
png_file <- file.path(out_dir, "fig_ar2_regions.png")

# The panel is fixed at 2:1 by coord_fixed, so the height follows from the width
# once the title, axis and legend bands are accounted for (about 0.75 in total).
fig_width  <- WIDTH_FRAC * TEXTWIDTH_IN
panel_w    <- fig_width - 0.52                     # minus y axis and margins
fig_height <- panel_w / 2 + 0.78

message(sprintf("Device size: %.2f x %.2f in", fig_width, fig_height))

if (!USE_TIKZ) {

  grDevices::pdf(pdf_file, width = fig_width, height = fig_height)
  print(fig)
  grDevices::dev.off()
  message("Preview PDF written to ", pdf_file)

} else {

  old_pkgs <- getOption("tikzLatexPackages")
  options(tikzLatexPackages = c(old_pkgs,
                                "\\usepackage{amsmath}\n",
                                "\\usepackage{amssymb}\n"))

  tikzDevice::tikz(
    file = tex_file, width = fig_width, height = fig_height,
    standAlone = TRUE, sanitize = FALSE, engine = "pdftex"
  )
  print(fig)
  grDevices::dev.off()

  options(tikzLatexPackages = old_pkgs)

  tex_src <- readLines(tex_file, warn = FALSE)
  n_paths <- sum(grepl("\\\\path|\\\\draw", tex_src))
  message("TikZ source: ", length(tex_src), " lines, ", n_paths,
          " path commands.")
  if (n_paths < 20L) {
    warning("Very few drawing commands: the panel is probably empty.")
  }

  pdflatex <- Sys.which("pdflatex")
  if (nzchar(pdflatex)) {
    old_wd <- setwd(out_dir)
    log_out <- system2(pdflatex,
                       args = c("-interaction=nonstopmode", basename(tex_file)),
                       stdout = TRUE, stderr = TRUE)
    setwd(old_wd)
    if (!file.exists(pdf_file)) {
      stop("LaTeX compilation failed:\n", paste(log_out, collapse = "\n"))
    }
    message("Figure PDF: ", pdf_file)
  } else {
    warning("pdflatex not found; the .tex was written but not compiled.")
  }
}

if (file.exists(pdf_file)) {
  if (requireNamespace("pdftools", quietly = TRUE)) {
    pdftools::pdf_convert(pdf = pdf_file, format = "png", pages = 1L,
                          dpi = 300, filenames = png_file, verbose = FALSE)
  } else {
    pdftoppm <- Sys.which("pdftoppm")
    if (nzchar(pdftoppm)) {
      system2(pdftoppm, args = c("-png", "-r", "300", "-singlefile", pdf_file,
                                 tools::file_path_sans_ext(png_file)))
    }
  }
}

if (file.exists(png_file)) message("Figure PNG: ", png_file)
