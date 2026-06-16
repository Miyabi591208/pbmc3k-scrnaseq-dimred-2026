#!/usr/bin/env Rscript

source("scripts/common.R")
ensure_dirs()

url <- "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
tar_path <- here("data", "pbmc3k_filtered_gene_bc_matrices.tar.gz")
extract_dir <- here("data")
matrix_dir <- here("data", "filtered_gene_bc_matrices", "hg19")

if (dir.exists(matrix_dir) && file.exists(file.path(matrix_dir, "matrix.mtx"))) {
  message("PBMC3k matrix already exists: ", matrix_dir)
  quit(save = "no", status = 0)
}

if (!file.exists(tar_path)) {
  message("Downloading PBMC3k from 10x Genomics...")
  ok <- tryCatch({
    download.file(url, destfile = tar_path, mode = "wb", quiet = FALSE)
    TRUE
  }, error = function(e) {
    message("Download failed: ", conditionMessage(e))
    FALSE
  })
  if (!ok) {
    stop(
      "Could not download PBMC3k. Check network access, then rerun this script.\nURL: ",
      url,
      call. = FALSE
    )
  }
}

message("Extracting: ", tar_path)
untar(tar_path, exdir = extract_dir)

if (!file.exists(file.path(matrix_dir, "matrix.mtx"))) {
  stop("Extraction finished, but matrix.mtx was not found at: ", matrix_dir, call. = FALSE)
}

message("PBMC3k is ready: ", matrix_dir)

