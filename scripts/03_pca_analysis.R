#!/usr/bin/env Rscript

source("scripts/common.R")
require_packages(c("Seurat", "ggplot2", "patchwork", "dplyr", "readr"))
ensure_dirs()
set_reproducible_seed()

pbmc <- load_rds_or_stop(
  here("results", "pbmc3k_preprocessed_pca.rds"),
  "Run scripts/02_seurat_preprocess.R first."
)

pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.5, random.seed = 20260616, verbose = FALSE)

p_pca <- DimPlot(pbmc, reduction = "pca", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("PCA scatter plot (colored by Seurat cluster)")
save_plot(p_pca, "pca_scatter_clusters.png", width = 7, height = 5)

p_elbow <- ElbowPlot(pbmc, ndims = 50) + ggtitle("Elbow plot: variance explained by PCs")
save_plot(p_elbow, "pca_elbow_plot.png", width = 7, height = 5)

p_load <- VizDimLoadings(pbmc, dims = 1:4, reduction = "pca", nfeatures = 12) +
  plot_annotation(title = "Top genes contributing to PC loadings")
save_plot(p_load, "pca_loading_top_genes_pc1_pc4.png", width = 12, height = 8)

p_heat <- DimHeatmap(pbmc, dims = 1:12, cells = 500, balanced = TRUE, fast = FALSE)
save_plot(p_heat, "pca_dimheatmap_pc1_pc12.png", width = 11, height = 12)

loadings <- Loadings(pbmc, reduction = "pca")[, 1:10]
write_csv(
  as.data.frame(loadings) %>% mutate(gene = rownames(loadings), .before = 1),
  here("results", "pca_loadings_pc1_pc10.csv")
)

saveRDS(pbmc, here("results", "pbmc3k_pca_clustered_dims10.rds"))
message("Saved clustered object: ", here("results", "pbmc3k_pca_clustered_dims10.rds"))

