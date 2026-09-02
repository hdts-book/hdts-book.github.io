library(ggplot2)
library(patchwork)
library(grid)

## =========================================================
## helpers
## =========================================================

make_sigma <- function(a = 0.85, b = 0.22, theta = pi / 4.5) {
  # major axis in positive direction (NE)
  R <- matrix(c(cos(theta), -sin(theta),
                sin(theta),  cos(theta)), 2, 2)
  R %*% diag(c(a^2, b^2), 2) %*% t(R)
}

mahal_sq <- function(x, mu, Sigma) {
  z <- x - mu
  as.numeric(t(z) %*% solve(Sigma, z))
}

ellipse_df <- function(mu, Sigma, level, n = 500) {
  ee <- eigen(Sigma, symmetric = TRUE)
  A  <- ee$vectors %*% diag(sqrt(ee$values), 2)
  tt <- seq(0, 2 * pi, length.out = n)
  xy <- t(mu + sqrt(level) * A %*% rbind(cos(tt), sin(tt)))
  data.frame(x = xy[, 1], y = xy[, 2])
}

diamond_df <- function(r = 1) {
  data.frame(
    x = c(-r, 0, r, 0),
    y = c(0, r, 0, -r)
  )
}

circle_df <- function(r = 1, n = 500) {
  tt <- seq(0, 2 * pi, length.out = n)
  data.frame(
    x = r * cos(tt),
    y = r * sin(tt)
  )
}

## =========================================================
## common settings
## =========================================================

xlim <- c(-1.45, 3.05)
ylim <- c(-1.15, 2.85)

Sigma <- make_sigma(a = 0.85, b = 0.22, theta = pi / 4.5)

# left: outer contour passes exactly through (0,1)
x_touch_lasso <- c(0, 1)
mu_lasso      <- c(0.95, 1.55)
outer_lasso   <- mahal_sq(x_touch_lasso, mu_lasso, Sigma)

# right: outer contour passes near (0,1), but NOT exactly there
x_touch_ridge <- c(0.18, sqrt(1 - 0.25^2))
mu_ridge      <- c(0.90, 1.60)
outer_ridge   <- mahal_sq(x_touch_ridge, mu_ridge, Sigma)

lev_lasso <- c(0.45, 0.72, 1.00) * outer_lasso
lev_ridge <- c(0.45, .6, 0.73) * outer_ridge

ell_lasso_1 <- ellipse_df(mu_lasso, Sigma, lev_lasso[1])
ell_lasso_2 <- ellipse_df(mu_lasso, Sigma, lev_lasso[2])
ell_lasso_3 <- ellipse_df(mu_lasso, Sigma, lev_lasso[3])   # outermost hits (0,1)

ell_ridge_1 <- ellipse_df(mu_ridge, Sigma, lev_ridge[1])
ell_ridge_2 <- ellipse_df(mu_ridge, Sigma, lev_ridge[2])
ell_ridge_3 <- ellipse_df(mu_ridge, Sigma, lev_ridge[3])   # outermost near (0,1), not exact

dia <- diamond_df(1)
cir <- circle_df(1)

theme_fig <- theme_void(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "#EBEBEB", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    plot.margin      = margin(4, 4, 4, 4)
  )

base_axes <- list(
  geom_segment(
    aes(x = xlim[1], y = 0, xend = xlim[2] - 0.12, yend = 0),
    arrow = arrow(length = unit(0.11, "inches")),
    linewidth = 0.35
  ),
  geom_segment(
    aes(x = 0, y = ylim[1], xend = 0, yend = ylim[2] - 0.12),
    arrow = arrow(length = unit(0.11, "inches")),
    linewidth = 0.35
  ),
  annotate("text", x = xlim[2] - 0.28, y = -0.10,
           label = "beta[1]", parse = TRUE, size = 4),
  annotate("text", x = -0.40, y = ylim[2] - 0.12,
           label = "beta[2]", parse = TRUE, size = 4)
)

## =========================================================
## left panel: lasso
## =========================================================

p_lasso <- ggplot() +
  coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) +
  geom_polygon(data = dia, aes(x, y), fill = "#00DFF0", color = NA, alpha = 0.9) +
  base_axes +
  geom_path(data = ell_lasso_1, aes(x, y), color = "#FF4D4D", linewidth = 0.45) +
  geom_path(data = ell_lasso_2, aes(x, y), color = "#FF4D4D", linewidth = 0.45) +
  geom_path(data = ell_lasso_3, aes(x, y), color = "#FF4D4D", linewidth = 0.45) +
  geom_point(aes(x = mu_lasso[1], y = mu_lasso[2]), size = 1.2) +
  annotate("text", x = mu_lasso[1] - 0.22, y = mu_lasso[2] - 0.16,
           label = "hat(beta)", parse = TRUE, size = 4) +
  theme_fig

## =========================================================
## right panel: ridge
## =========================================================

p_ridge <- ggplot() +
  coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) +
  geom_polygon(data = cir, aes(x, y), fill = "#00DFF0", color = NA, alpha = 0.9) +
  base_axes +
  geom_path(data = ell_ridge_1, aes(x, y), color = "#FF4D4D", linewidth = 0.45) +
  geom_path(data = ell_ridge_2, aes(x, y), color = "#FF4D4D", linewidth = 0.45) +
  geom_path(data = ell_ridge_3, aes(x, y), color = "#FF4D4D", linewidth = 0.45) +
  geom_point(aes(x = mu_ridge[1], y = mu_ridge[2]), size = 1.2) +
  annotate("text", x = mu_ridge[1] - 0.22, y = mu_ridge[2] - 0.16,
           label = "hat(beta)", parse = TRUE, size = 4) +
  theme_fig

p_all <- p_lasso + p_ridge + plot_layout(ncol = 2)

print(p_all)

ggsave("lasso_ridge_corrected.png",
       p_all, width = 6.4, height = 3.0, dpi = 300)

ggsave("lasso_ridge_corrected.pdf",
       p_all, width = 6.4, height = 3.0,
       device = cairo_pdf)
