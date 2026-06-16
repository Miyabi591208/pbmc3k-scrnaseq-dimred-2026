suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
})

project_root <- function() {
  normalizePath(file.path(dirname(sys.frame(1)$ofile %||% getwd()), ".."), mustWork = FALSE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

here <- function(...) {
  root <- normalizePath(file.path(getwd()), mustWork = FALSE)
  if (basename(root) == "scripts") root <- normalizePath(file.path(root, ".."), mustWork = FALSE)
  file.path(root, ...)
}

ensure_dirs <- function() {
  dirs <- c("data", "results", "figures", "gifs")
  for (d in dirs) dir.create(here(d), showWarnings = FALSE, recursive = TRUE)
}

set_reproducible_seed <- function(seed = 20260616) {
  set.seed(seed)
  options(Seurat.object.assay.version = "v5")
}

save_plot <- function(plot, filename, width = 7, height = 5, dpi = 160) {
  ensure_dirs()
  path <- here("figures", filename)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  message("saved: ", path)
  invisible(path)
}

require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      "\nRun scripts/00_setup.R first, or install them manually.",
      call. = FALSE
    )
  }
}

load_rds_or_stop <- function(path, hint) {
  if (!file.exists(path)) stop("Missing file: ", path, "\n", hint, call. = FALSE)
  readRDS(path)
}

