# seurat_clustering: A quick end-to-end clustering of a Seurat object

Wrapper to perform PCA, neighbor-finding, UMAP calculations and
clustering

## Usage

``` r
seurat_clustering(so_in, npcs_in, resolution = 0.8, algorithm = 3)
```

## Arguments

- so_in:

  A Seurat single cell object

- npcs_in:

  Number of principal components to be used to calculate neighbors, UMAP
  projection, and clustering

- resolution:

  A resolution value for clustering. Higher values will produce more
  clusters (Default 0.8)

- algorithm:

  A value of 1 (Louvain), 2 (Leiden), or 3 (Slow Local Moving (SLM)) to
  select the clustering algorithm (Default 3)

## Value

The updated Seurat single cell object with recalculated PCA, neighbors,
UMAP and clusters
