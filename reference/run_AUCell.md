# run_AUCell: A wrapper script for running AUCell for per cell gene set enrichment analysis

Runs AUCell on basic inputs

## Usage

``` r
run_AUCell(so, gene_sets)
```

## Arguments

- so:

  A Seurat object

- gene_sets:

  A list of gene sets, such as those obtained via GSEABase

## Value

A list object containing the AUCell scores for each gene set and the
assignment of the cells as to the likelihood of the cells overexpressing
the gene set. These values can be added to the original Seurat object ex
post facto.

## Details

AUCell evaluates the likelihood of a gene set being enriched in an
individual cells using the Area Under the Curve (AUC). The wrapper
script uses a Seurat object and a list of gene sets to calculate the AUC
for each gene set and also provides the assignment likelihood, based on
the distribution of the AUCs across the dataset.
