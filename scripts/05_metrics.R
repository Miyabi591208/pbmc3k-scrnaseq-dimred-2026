#!/usr/bin/env Rscript

source("scripts/common.R")
require_packages(c("Seurat", "FNN", "mclust", "aricode", "dplyr", "readr"))
ensure_dirs()
set_reproducible_seed()

knn_preservation <- function(high, low, k = 15) {
  high_nn <- FNN::get.knn(high, k = k)$nn.index
  low_nn <- FNN::get.knn(low, k = k)$nn.index
  mean(vapply(seq_len(nrow(high_nn)), function(i) {
    length(intersect(high_nn[i, ], low_nn[i, ])) / k
  }, numeric(1)))
}

trustworthiness_score <- function(high, low, k = 15) {
  n <- nrow(high)
  if (n <= 2 * k + 1) stop("n must be greater than 2k + 1 for trustworthiness.")
  high_dist <- as.matrix(dist(high))
  low_nn <- FNN::get.knn(low, k = k)$nn.index
  high_rank <- t(apply(high_dist, 1, rank, ties.method = "average"))
  penalty <- 0
  for (i in seq_len(n)) {
    for (j in low_nn[i, ]) {
      r <- high_rank[i, j]
      if (r > k) penalty <- penalty + (r - k)
    }
  }
  1 - (2 / (n * k * (2 * n - 3 * k - 1))) * penalty
}

marker_sets <- list(
  T_NK = c("CD3D", "CD3E", "IL7R", "NKG7", "GNLY"),
  B = c("MS4A1", "CD79A"),
  Mono = c("LYZ", "S100A8", "S100A9", "FCGR3A", "MS4A7"),
  DC = c("FCER1A", "CST3"),
  Platelet = c("PPBP", "PF4")
)

files <- list.files(here("results"), pattern = "^pbmc3k_dims[0-9]+_umap_tsne_clustered\\.rds$", full.names = TRUE)
if (length(files) == 0) {
  stop("No sweep RDS files found. Run scripts/04_tsne_umap_parameter_sweep.R first.", call. = FALSE)
}

rows <- list()
for (f in files) {
  obj <- readRDS(f)
  dims <- as.integer(sub(".*dims([0-9]+)_.*", "\\1", basename(f)))
  high <- Embeddings(obj, "pca")[, seq_len(dims), drop = FALSE]
  umap <- Embeddings(obj, paste0("umap_dims", dims))
  tsne <- Embeddings(obj, paste0("tsne_dims", dims))
  clusters <- obj$seurat_clusters
  ref <- readRDS(here("results", "pbmc3k_dims20_umap_tsne_clustered.rds"))$seurat_clusters

  rows[[paste0("umap_", dims)]] <- tibble(
    method = "UMAP",
    dims = dims,
    trustworthiness_k15 = trustworthiness_score(high, umap, k = 15),
    knn_preservation_k15 = knn_preservation(high, umap, k = 15),
    cluster_count = length(unique(clusters)),
    ari_vs_dims20_cluster = mclust::adjustedRandIndex(as.integer(clusters), as.integer(ref)),
    nmi_vs_dims20_cluster = aricode::NMI(as.integer(clusters), as.integer(ref))
  )
  rows[[paste0("tsne_", dims)]] <- tibble(
    method = "t-SNE",
    dims = dims,
    trustworthiness_k15 = trustworthiness_score(high, tsne, k = 15),
    knn_preservation_k15 = knn_preservation(high, tsne, k = 15),
    cluster_count = length(unique(clusters)),
    ari_vs_dims20_cluster = mclust::adjustedRandIndex(as.integer(clusters), as.integer(ref)),
    nmi_vs_dims20_cluster = aricode::NMI(as.integer(clusters), as.integer(ref))
  )
}

metrics <- bind_rows(rows)
write_csv(metrics, here("results", "embedding_metrics.csv"))

obj <- readRDS(here("results", "pbmc3k_dims20_umap_tsne_clustered.rds"))
avg <- AverageExpression(obj, features = unique(unlist(marker_sets)), assays = "RNA", slot = "data")$RNA
marker_summary <- lapply(names(marker_sets), function(cell_type) {
  genes <- intersect(marker_sets[[cell_type]], rownames(avg))
  if (length(genes) == 0) return(NULL)
  data.frame(cell_type = cell_type, gene = genes, avg[genes, , drop = FALSE], check.names = FALSE)
}) %>% bind_rows()
write_csv(marker_summary, here("results", "marker_average_expression_by_cluster.csv"))

message("Saved metrics: ", here("results", "embedding_metrics.csv"))
