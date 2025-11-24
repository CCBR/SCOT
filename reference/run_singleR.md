# run_singleR: Runs SingleR directly on the Seurat object

Wrapper script for the SingleR function

## Usage

``` r
run_singleR(so, ref_file, label)
```

## Arguments

- so:

  A Seurat single cell object

- ref_file:

  A reference file in a SingleCellExperiment format, such as those
  obtained from SingleR/CellDex package

- label:

  The label identity to be used. Must be a column header in the metadata
  of `ref_file`

## Value

A vector of pruned cell type labels
