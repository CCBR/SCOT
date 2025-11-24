# add_VDJ_metadata: Appends V(D)J alignments from CellRanger to an existing Seurat object

Wrapper for djvdj to assign V(D)J alignments to cells in a Seurat
object. Can be used for TCRa/TCRb, TCRg/TCRd, and BCR

## Usage

``` r
add_VDJ_metadata(so, vdj_file)
```

## Arguments

- so:

  A Seurat single cell RNA object

- vdj_file:

  A path to a V(D)J annotation directory from CellRanger

## Value

A Seurat object with V(D)J metadata
