#!/usr/bin/env Rscript

source("scripts/common.R")
require_packages(c("Seurat", "ggplot2", "patchwork", "dplyr", "readr"))
ensure_dirs()
set_reproducible_seed()

pbmc_base <- load_rds_or_stop(
  here("results", "pbmc3k_preprocessed_pca.rds"),
  "Run scripts/02_seurat_preprocess.R first."
)

dims_list <- c(5, 10, 20, 30, 50)
perplexities <- c(5, 10, 30, 50, 100)
neighbors <- c(5, 15, 30, 50)
min_dists <- c(0.01, 0.1, 0.3, 0.8)

timings <- list()
cluster_summary <- list()

for (d in dims_list) {
  dims_use <- seq_len(d)
  message("PC sweep dims=1:", d)
  obj <- pbmc_base
  t0 <- Sys.time()
  obj <- FindNeighbors(obj, dims = dims_use, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, random.seed = 20260616, verbose = FALSE)
  cluster_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  t1 <- Sys.time()
  obj <- RunUMAP(
    obj,
    dims = dims_use,
    reduction.name = paste0("umap_dims", d),
    reduction.key = paste0("UMAPD", d, "_"),
    n.neighbors = 30,
    min.dist = 0.3,
    metric = "cosine",
    seed.use = 20260616,
    verbose = FALSE
  )
  umap_time <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

  t2 <- Sys.time()
  obj <- RunTSNE(
    obj,
    dims = dims_use,
    reduction.name = paste0("tsne_dims", d),
    reduction.key = paste0("TSNED", d, "_"),
    perplexity = 30,
    seed.use = 20260616,
    check_duplicates = FALSE,
    verbose = FALSE
  )
  tsne_time <- as.numeric(difftime(Sys.time(), t2, units = "secs"))

  p_umap <- DimPlot(obj, reduction = paste0("umap_dims", d), group.by = "seurat_clusters", label = TRUE) +
    ggtitle(paste0("UMAP: dims=1:", d, ", n_neighbors=30, min_dist=0.3"))
  p_tsne <- DimPlot(obj, reduction = paste0("tsne_dims", d), group.by = "seurat_clusters", label = TRUE) +
    ggtitle(paste0("t-SNE: dims=1:", d, ", perplexity=30"))
  save_plot(p_umap, sprintf("umap_dims%d_neighbors30_mindist0.3.png", d), width = 7, height = 5)
  save_plot(p_tsne, sprintf("tsne_dims%d_perplexity30.png", d), width = 7, height = 5)

  timings[[paste0("dims", d)]] <- tibble(
    analysis = "pc_sweep",
    dims = d,
    perplexity = NA_real_,
    n_neighbors = 30,
    min_dist = 0.3,
    metric = "cosine",
    cluster_seconds = cluster_time,
    umap_seconds = umap_time,
    tsne_seconds = tsne_time
  )
  cluster_summary[[paste0("dims", d)]] <- tibble(
    analysis = "pc_sweep",
    dims = d,
    cluster_count = length(unique(obj$seurat_clusters))
  )
  saveRDS(obj, here("results", sprintf("pbmc3k_dims%d_umap_tsne_clustered.rds", d)))
}

pbmc_ref <- pbmc_base
pbmc_ref <- FindNeighbors(pbmc_ref, dims = 1:20, verbose = FALSE)
pbmc_ref <- FindClusters(pbmc_ref, resolution = 0.5, random.seed = 20260616, verbose = FALSE)

tsne_plots <- list()
for (p in perplexities) {
  message("t-SNE perplexity=", p)
  obj <- RunTSNE(
    pbmc_ref,
    dims = 1:20,
    reduction.name = paste0("tsne_p", p),
    reduction.key = paste0("TSNEP", p, "_"),
    perplexity = p,
    seed.use = 20260616,
    check_duplicates = FALSE,
    verbose = FALSE
  )
  plot <- DimPlot(obj, reduction = paste0("tsne_p", p), group.by = "seurat_clusters", label = FALSE) +
    ggtitle(paste0("perplexity=", p))
  tsne_plots[[paste0("p", p)]] <- plot
  save_plot(plot, sprintf("tsne_dims20_perplexity%d.png", p), width = 6, height = 5)
}
save_plot(wrap_plots(tsne_plots, ncol = 3), "tsne_perplexity_grid_dims20.png", width = 14, height = 9)

neighbor_plots <- list()
for (nn in neighbors) {
  message("UMAP n_neighbors=", nn)
  obj <- RunUMAP(
    pbmc_ref,
    dims = 1:20,
    reduction.name = paste0("umap_nn", nn),
    reduction.key = paste0("UMAPNN", nn, "_"),
    n.neighbors = nn,
    min.dist = 0.3,
    metric = "cosine",
    seed.use = 20260616,
    verbose = FALSE
  )
  plot <- DimPlot(obj, reduction = paste0("umap_nn", nn), group.by = "seurat_clusters", label = FALSE) +
    ggtitle(paste0("n_neighbors=", nn, ", min_dist=0.3"))
  neighbor_plots[[paste0("nn", nn)]] <- plot
  save_plot(plot, sprintf("umap_dims20_neighbors%d_mindist0.3.png", nn), width = 6, height = 5)
}
save_plot(wrap_plots(neighbor_plots, ncol = 2), "umap_n_neighbors_grid_dims20.png", width = 11, height = 9)

mindist_plots <- list()
for (md in min_dists) {
  md_label <- gsub("\\.", "p", as.character(md))
  message("UMAP min_dist=", md)
  obj <- RunUMAP(
    pbmc_ref,
    dims = 1:20,
    reduction.name = paste0("umap_md", md_label),
    reduction.key = paste0("UMAPMD", md_label, "_"),
    n.neighbors = 30,
    min.dist = md,
    metric = "cosine",
    seed.use = 20260616,
    verbose = FALSE
  )
  plot <- DimPlot(obj, reduction = paste0("umap_md", md_label), group.by = "seurat_clusters", label = FALSE) +
    ggtitle(paste0("n_neighbors=30, min_dist=", md))
  mindist_plots[[paste0("md", md_label)]] <- plot
  save_plot(plot, sprintf("umap_dims20_neighbors30_mindist%s.png", md_label), width = 6, height = 5)
}
save_plot(wrap_plots(mindist_plots, ncol = 2), "umap_min_dist_grid_dims20.png", width = 11, height = 9)

grid_plots <- list()
for (nn in neighbors) {
  for (md in min_dists) {
    md_label <- gsub("\\.", "p", as.character(md))
    message("UMAP grid nn=", nn, " min_dist=", md)
    obj <- RunUMAP(
      pbmc_ref,
      dims = 1:20,
      reduction.name = paste0("umap_nn", nn, "_md", md_label),
      reduction.key = paste0("UMAPG", nn, md_label, "_"),
      n.neighbors = nn,
      min.dist = md,
      metric = "cosine",
      seed.use = 20260616,
      verbose = FALSE
    )
    plot <- DimPlot(obj, reduction = paste0("umap_nn", nn, "_md", md_label), group.by = "seurat_clusters", label = FALSE) +
      ggtitle(paste0("nn=", nn, ", md=", md)) +
      NoLegend()
    grid_plots[[paste0("nn", nn, "_md", md_label)]] <- plot
  }
}
save_plot(wrap_plots(grid_plots, ncol = 4), "umap_neighbors_min_dist_grid_dims20.png", width = 16, height = 14)

write_csv(bind_rows(timings), here("results", "runtime_summary.csv"))
write_csv(bind_rows(cluster_summary), here("results", "cluster_summary.csv"))
message("Parameter sweep complete.")
