#!/usr/bin/env Rscript

source("scripts/common.R")
require_packages(c("Seurat", "ggplot2", "patchwork", "dplyr", "readr"))
ensure_dirs()
set_reproducible_seed()

matrix_dir <- here("data", "filtered_gene_bc_matrices", "hg19")
if (!dir.exists(matrix_dir)) {
  stop("PBMC3k data not found. Run scripts/01_download_pbmc3k.R first.", call. = FALSE)
}

message("Reading 10x matrix: ", matrix_dir)
pbmc.data <- Read10X(data.dir = matrix_dir)

pbmc <- CreateSeuratObject(
  counts = pbmc.data,
  project = "pbmc3k",
  min.cells = 3,
  min.features = 200
)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

qc_before <- pbmc@meta.data %>%
  mutate(cell = rownames(pbmc@meta.data))
write_csv(qc_before, here("results", "qc_metrics_before_filter.csv"))

p_qc_before <- VlnPlot(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.05
) + plot_annotation(title = "PBMC3k QC metrics before filtering")
save_plot(p_qc_before, "qc_violin_before_filter.png", width = 10, height = 4)

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)
p_var <- VariableFeaturePlot(pbmc)
p_var_label <- LabelPoints(plot = p_var, points = top10, repel = TRUE)
save_plot(p_var + p_var_label, "variable_features_top10.png", width = 11, height = 5)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, verbose = FALSE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 50, verbose = FALSE)

saveRDS(pbmc, here("results", "pbmc3k_preprocessed_pca.rds"))
message("Saved preprocessed object: ", here("results", "pbmc3k_preprocessed_pca.rds"))

