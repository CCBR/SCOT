# split_featurePlot: Splits FeaturePlot by groups

Adapts the FeaturePlot function from Seurat to split by groups and
arrange by user specification

## Usage

``` r
split_featurePlot(
  so,
  features,
  split_ident,
  label = FALSE,
  ncol = NA,
  nrow = NA,
  min.cutoff = NA,
  max.cutoff = NA,
  plot_image = FALSE,
  return_list = FALSE,
  slot = "scale.data",
  order = FALSE,
  reduction = NULL
)
```

## Arguments

- so:

  A Seurat object

- features:

  A list of features (e.g. genes) to plot on a projection

- split_ident:

  A metadata identity for splitting the image

- label:

  Boolean for whether to label the plot with the current Idents value of
  the Seurat object

- ncol:

  Integer value for number of columns in plot

- nrow:

  Integer value for number of rows in plot

- min.cutoff:

  Vector of minimum cutoff values for each feature, may specify quantile
  in the form of 'q##' where '##' is the quantile (eg, 1, 10)

- max.cutoff:

  Vector of maximum cutoff values for each feature, may specify quantile
  in the form of 'q##' where '##' is the quantile (eg, 1, 10)

- plot_image:

  Boolean for whether to export plots to current visualization window

- return_list:

  Boolean for whether to export list of plots

- slot:

  Data slot from which to extract values

- order:

  Boolean for whether to show highest feature values at the front of the
  image

- reduction:

  Dimensionality reduction to use

## Value

A list of ggplot2 plots
