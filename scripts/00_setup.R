#!/usr/bin/env Rscript

required_cran <- c(
  "Seurat", "SeuratObject", "ggplot2", "patchwork", "dplyr", "readr",
  "future", "FNN", "mclust", "aricode", "magick"
)

installed <- rownames(installed.packages())
missing <- setdiff(required_cran, installed)

if (length(missing) == 0) {
  message("All required R packages are already installed.")
} else {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

message("R version: ", R.version.string)
message("Seurat version: ", as.character(packageVersion("Seurat")))

