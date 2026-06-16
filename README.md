# PBMC3k dimensionality reduction article

This project contains a reproducible analysis for the article:

`【2026年版】scRNA-seq の次元削減を数式から理解する：PCA・t-SNE・UMAP をPBMC3kで実際に動かして比較`

The analysis uses the public 10x Genomics PBMC3k dataset and compares PCA, t-SNE, and UMAP with controlled random seeds.

## Data

- Dataset: 10x Genomics, `3k PBMCs from a Healthy Donor`
- Source: https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0
- Cells: 2,700 detected cells
- License: CC BY 4.0, as stated on the 10x Genomics dataset page
- Tutorial compatibility: same PBMC3k dataset family used in the Seurat guided clustering tutorial

The data archive is not intended to be committed. It is downloaded by `scripts/01_download_pbmc3k.R`.

## Environment

Tested locally with:

- R 4.4.2
- Seurat 5.4.0
- SeuratObject 5.0.2

Install R packages:

```bash
Rscript scripts/00_setup.R
```

Required R packages:

- Seurat
- SeuratObject
- ggplot2
- patchwork
- dplyr
- readr
- future
- FNN
- mclust
- aricode
- magick

Python is optional. The included GIF script uses R + magick so that the core workflow remains R-based.

If you want to create alternative Python animations, use a virtual environment with:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install numpy pandas matplotlib pillow scikit-learn umap-learn
```

## Run Order

Run commands from the project root.

```bash
Rscript scripts/00_setup.R
Rscript scripts/01_download_pbmc3k.R
Rscript scripts/02_seurat_preprocess.R
Rscript scripts/03_pca_analysis.R
Rscript scripts/04_tsne_umap_parameter_sweep.R
Rscript scripts/05_metrics.R
Rscript scripts/06_make_gifs.R
```

## Outputs

Figures are written to `figures/`.

Expected figure files include:

- `qc_violin_before_filter.png`
- `variable_features_top10.png`
- `pca_scatter_clusters.png`
- `pca_elbow_plot.png`
- `pca_loading_top_genes_pc1_pc4.png`
- `pca_dimheatmap_pc1_pc12.png`
- `umap_dims5_neighbors30_mindist0.3.png`
- `umap_dims10_neighbors30_mindist0.3.png`
- `umap_dims20_neighbors30_mindist0.3.png`
- `umap_dims30_neighbors30_mindist0.3.png`
- `umap_dims50_neighbors30_mindist0.3.png`
- `tsne_dims5_perplexity30.png`
- `tsne_dims10_perplexity30.png`
- `tsne_dims20_perplexity30.png`
- `tsne_dims30_perplexity30.png`
- `tsne_dims50_perplexity30.png`
- `tsne_perplexity_grid_dims20.png`
- `umap_n_neighbors_grid_dims20.png`
- `umap_min_dist_grid_dims20.png`
- `umap_neighbors_min_dist_grid_dims20.png`

GIFs are written to `gifs/`.

Expected GIF files include:

- `pca_axis_rotation_toy.gif`
- `tsne_iterative_layout_toy.gif`
- `umap_neighbors_intuition_toy.gif`
- `pbmc3k_umap_min_dist_sweep.gif`

Tables and intermediate RDS files are written to `results/`.

Important result files include:

- `pbmc3k_preprocessed_pca.rds`
- `pbmc3k_pca_clustered_dims10.rds`
- `runtime_summary.csv`
- `cluster_summary.csv`
- `embedding_metrics.csv`
- `marker_average_expression_by_cluster.csv`

## Reproducibility Notes

- Random seed is fixed to `20260616`.
- Seurat UMAP and t-SNE are stochastic; exact coordinates can differ across Seurat, uwot, Rtsne, BLAS, and operating system versions.
- The biological interpretation should not rely on a single embedding.
- Cluster IDs are convenience labels, not ground-truth cell-type labels.
- PBMC marker checks are included as sanity checks, not as a complete annotation workflow.

## Runtime Estimate

On a recent laptop, the preprocessing and PCA steps typically finish in minutes. The parameter sweep is slower because it runs multiple UMAP and t-SNE configurations; expect tens of minutes depending on CPU and package versions.

## Citation Targets

Please cite:

- 10x Genomics PBMC3k dataset
- Seurat guided clustering tutorial
- van der Maaten and Hinton, 2008
- Wattenberg, Viégas and Johnson, 2016
- Kobak and Berens, 2019
- McInnes, Healy and Melville, UMAP
- scikit-learn trustworthiness documentation if using the metric explanation

