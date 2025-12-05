# make_bubble_plot: Flexible dotplot

Creates a bubble plot that shows average expression percent expression
for genes in single cell categories

## Usage

``` r
make_bubble_plot(
  so,
  features,
  palette = "RdBu",
  assay = "SCT",
  scale = FALSE,
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

- scale:

  Boolean to scale average expression across identities

- ident:

  The categorical identity to classify groups of cells

## Value

A ggplot2 figure
