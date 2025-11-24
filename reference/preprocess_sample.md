# preprocess_sample: Performs an initial preprocess on a Seurat object

Conducts preliminary normalization/scaling, cell cycle scoring, and
clustering on a Seurat single cell object

## Usage

``` r
preprocess_sample(so_in, species, npcs_in)
```

## Arguments

- so_in:

  A Seurat single cell object with RNA counts

- species:

  Species of the sample (hg38, hg19, or mm10)

- npcs_in:

  Number of principal components to use for dimensionality reduction

## Value

A Seurat single cell object with normalized read counts

## Details

The function performs the following processes on the data:assay

- Cell cycle scoring

- Normalization and scaling of RNA counts with log normalization

- Normalization and scaling of RNA counts with SCTransform

- Principal components analysis, neighbor identification, UMAP
  dimensionality reduction and preliminary clustering (through
  `seurat_clustering` function)
