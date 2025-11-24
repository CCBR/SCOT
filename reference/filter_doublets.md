# filter_doublets: Identifies and removes doublets or multiplets in scRNA samples

Wrapper for DoubletFinder and scDblFinder for identifying doublets

## Usage

``` r
filter_doublets(so_in, doublet_finder_method = "DoubletFinder")
```

## Arguments

- so_in:

  A provided Seurat scRNA object

- doublet_finder_method:

  Character string indicating the use of "DoubletFinder" (default),
  "scDblFinder", "consensus" removal of doublets, or the "union" of
  total doublets identified in both algorithms

## Value

A subsetted Seurat object with doublets removed

## Details

Runs the latest version of DoubletFinder and/or scDblFinder and returns
a filtered scRNA object with doublets removed
