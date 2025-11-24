# run_batch_correction: Performs batch correction within Seurat v5 framework

Performs batch correction using Seurat v5, with additional
dimensionality reduction (UMAP), clustering, and cluster-based cell type
annotation.

## Usage

``` r
run_batch_correction(
  so_in,
  npcs,
  species,
  resolution_list,
  method_in,
  reduction_in = c(0.2, 0.4, 0.5, 0.8, 1),
  vars_to_regress = NULL,
  conda_env = ""
)
```

## Arguments

- so_in:

  A merged Seurat object containing all samples prior to batch
  correction

- npcs:

  The number of principal components to use for dimensionality reduction
  and neighbor identification

- species:

  "hg19", "hg38", or "mm10". Used to determine databases used for cell
  type annotation

- resolution_list:

  A vector of resolutions ranging from 0 to 2.0 for clustering. Smaller
  resolutions produce fewer clusters (Default)

- method_in:

  A character string to indicate which batch correction method to use

- reduction_in:

  A character string for naming the PCA produced

- vars_to_regress:

  A character vector for variables to regress out when running
  SCTransform normalization and batch correction

- conda_env:

  A character string for indicating which conda environment contains the
  necessary packages

## Value

A batch corrected Seurat object with updated UMAP projections and
clusters

## Details

Utilizes the Seurat v5 framework to streamline batch correction using
one of the following options:

- SCVI

- LIGER

- Canonical Correlation Analysis (CCA, or Seurat "integrate")

- Reciprocal PCA (RPCA)

- Harmony Follows up with re-running dimensionality reduction with PCA,
  neighbor finding, and UMAP projection. Also runs clustering with the
  slow local moving algorithm at a series of resolutions selected by the
  user. Final step conducts cell type annotation with SingleR using the
  average gene expression vector of each cluster
