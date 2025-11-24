# make a flexible bubble plot with Seurat::DotPlot

Creates a bubble plot that shows average expression percent expression
for genes in single cell categories

## Usage

``` r
make_bubble_plot(
  so,
  features,
  palette = "RdBu",
  assay = "SCT",
  ident = "seurat_clusters"
)
```

## Arguments

- so:

  A Seurat single cell RNA object

- features:

  A character vector of genes to plot

- palette:

  A color palette from ggplot2

- assay:

  The counts assay to use for determining expression

- ident:

  The categorical identity to classify groups of cells

## Value

A ggplot2 figure
