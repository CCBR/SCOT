# run_singleR_cluster: Runs SingleR on average expression of clusters

Wrapper function to run SingleR cell type annotation on the current
default identities of a Seurat scRNA object, e.g. clusters

## Usage

``` r
run_singleR_cluster(so_in, ref_file, label)
```

## Arguments

- so_in:

  A Seurat single cell object

- ref_file:

  A SingleR compatible single cell reference object
  (SingleCellExperiment object)

- label:

  A cell type label from the metadata column headers

## Value

A character vector of matched cell type annotations based on clusters
