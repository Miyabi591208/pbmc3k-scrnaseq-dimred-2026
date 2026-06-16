#!/usr/bin/env Rscript

source("scripts/common.R")
require_packages(c("ggplot2", "magick"))
ensure_dirs()
set_reproducible_seed()

make_gif <- function(frame_files, out_file, delay = 12) {
  img <- magick::image_read(frame_files)
  info <- magick::image_info(img)
  canvas <- sprintf("%dx%d", max(info$width), max(info$height))
  img <- magick::image_extent(img, geometry = canvas, gravity = "center", color = "white")
  anim <- magick::image_animate(img, delay = delay, loop = 0, dispose = "background", optimize = FALSE)
  magick::image_write(anim, out_file)
  message("saved: ", out_file)
}

tmp_dir <- here("results", "gif_frames")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

toy <- data.frame(
  x = c(rnorm(120, -1.2, 0.45), rnorm(120, 1.2, 0.45)),
  y = c(rnorm(120, -0.2, 0.25), rnorm(120, 0.4, 0.25)),
  group = rep(c("A", "B"), each = 120)
)
theta <- pi / 5
rot <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
toy[, c("x", "y")] <- as.matrix(toy[, c("x", "y")]) %*% rot
pc <- prcomp(toy[, c("x", "y")], center = TRUE, scale. = TRUE)$rotation[, 1]

frames <- character()
for (i in seq_len(24)) {
  alpha <- (i - 1) / 23
  angle <- atan2(pc[2], pc[1]) * alpha
  axis_end <- data.frame(x = c(0, 2.4 * cos(angle)), y = c(0, 2.4 * sin(angle)))
  p <- ggplot(toy, aes(x, y, color = group)) +
    geom_point(alpha = 0.65, size = 1.8) +
    geom_path(data = axis_end, aes(x, y), inherit.aes = FALSE, linewidth = 1.4, color = "black") +
    coord_equal(xlim = c(-3, 3), ylim = c(-3, 3)) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none") +
    labs(title = "PCA: variance-maximizing axis rotates into place", x = "gene axis 1", y = "gene axis 2")
  frame <- file.path(tmp_dir, sprintf("pca_axis_%02d.png", i))
  ggsave(frame, p, width = 5, height = 5, dpi = 110, bg = "white")
  frames <- c(frames, frame)
}
make_gif(frames, here("gifs", "pca_axis_rotation_toy.gif"), delay = 10)

base <- toy
frames <- character()
for (i in seq_len(24)) {
  alpha <- (i - 1) / 23
  out <- base
  out$x <- (1 - alpha) * rnorm(nrow(base), 0, 1.3) + alpha * base$x
  out$y <- (1 - alpha) * rnorm(nrow(base), 0, 1.3) + alpha * base$y
  p <- ggplot(out, aes(x, y, color = group)) +
    geom_point(alpha = 0.75, size = 1.8) +
    coord_equal(xlim = c(-3.2, 3.2), ylim = c(-3.2, 3.2)) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none") +
    labs(title = sprintf("t-SNE intuition: iterative layout, step %02d", i), x = "t-SNE 1", y = "t-SNE 2")
  frame <- file.path(tmp_dir, sprintf("tsne_iter_%02d.png", i))
  ggsave(frame, p, width = 5, height = 5, dpi = 110, bg = "white")
  frames <- c(frames, frame)
}
make_gif(frames, here("gifs", "tsne_iterative_layout_toy.gif"), delay = 10)

frames <- character()
for (i in seq_len(24)) {
  alpha <- (i - 1) / 23
  out <- base
  global <- scale(base[, c("x", "y")])[, 1:2]
  local <- base[, c("x", "y")]
  out$x <- (1 - alpha) * local$x + alpha * as.numeric(global[, 1]) * 1.2
  out$y <- (1 - alpha) * local$y + alpha * as.numeric(global[, 2]) * 1.2
  p <- ggplot(out, aes(x, y, color = group)) +
    geom_point(alpha = 0.75, size = 1.8) +
    coord_equal(xlim = c(-3.2, 3.2), ylim = c(-3.2, 3.2)) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none") +
    labs(title = sprintf("UMAP intuition: local to broader neighborhoods, frame %02d", i), x = "UMAP 1", y = "UMAP 2")
  frame <- file.path(tmp_dir, sprintf("umap_neighbors_%02d.png", i))
  ggsave(frame, p, width = 5, height = 5, dpi = 110, bg = "white")
  frames <- c(frames, frame)
}
make_gif(frames, here("gifs", "umap_neighbors_intuition_toy.gif"), delay = 10)

grid_file <- here("figures", "umap_min_dist_grid_dims20.png")
if (file.exists(grid_file)) {
  img <- magick::image_read(grid_file)
  info <- magick::image_info(img)
  cell_width <- floor(info$width[1] / 2)
  cell_height <- floor(info$height[1] / 2)
  crops <- c(
    sprintf("%dx%d+0+0", cell_width, cell_height),
    sprintf("%dx%d+%d+0", cell_width, cell_height, cell_width),
    sprintf("%dx%d+0+%d", cell_width, cell_height, cell_height),
    sprintf("%dx%d+%d+%d", cell_width, cell_height, cell_width, cell_height)
  )
  frames <- character()
  for (i in seq_along(crops)) {
    frame <- file.path(tmp_dir, sprintf("pbmc_mindist_%02d.png", i))
    cropped <- magick::image_crop(img, crops[i])
    cropped <- magick::image_extent(cropped, geometry = sprintf("%dx%d", cell_width, cell_height), gravity = "center", color = "white")
    magick::image_write(cropped, frame)
    frames <- c(frames, frame)
  }
  make_gif(frames, here("gifs", "pbmc3k_umap_min_dist_sweep.gif"), delay = 70)
} else {
  message("Skipping PBMC3k min_dist GIF because figure is missing: ", grid_file)
}

message("GIF generation complete.")
