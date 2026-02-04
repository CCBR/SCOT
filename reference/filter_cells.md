# filter_cells: A wrapper for multiple scRNA filtering methods

Wrapper function to filter cells from a Seurat object using method of
choice. Returns a Seurat object with filtering results

## Usage

``` r
filter_cells(
  so,
  method = "miQC",
  mads = 3,
  nCount_RNA_max = 50000,
  nCount_RNA_min = 1000,
  nFeature_RNA_max = 5000,
  nFeature_RNA_min = 200,
  percent_mt_max = 10,
  percent_mt_min = 0
)
```

## Arguments

- so:

  Seurat object to filter

- method:

  Filtering method to use: "miQC" (default), "manual", or "mads"

- mads:

  Number of median absolute deviations to use for "mads" method (default
  3)

- nCount_RNA_max:

  Maximum nCount_RNA threshold for "manual" method (default 500000)

- nCount_RNA_min:

  Minimum nCount_RNA threshold for "manual" method (default 1000)

- nFeature_RNA_max:

  Maximum nFeature_RNA threshold for "manual" method (default 5000)

- nFeature_RNA_min:

  Minimum nFeature_RNA threshold for "manual" method (default 200)

- percent_mt_max:

  Maximum percent_mt threshold for "manual" method (default 10)

- percent_mt_min:

  Minimum percent_mt threshold for "manual" method (default 0)

## Value

Returns a Seurat object with an additional metadata column "keep"
indicating whether each cell is to be kept or discarded

## Details

This function filters cells from a Seurat object using one of three
methods:

- "miQC": Uses the miQC package to filter cells based on a probabilistic
  model based on mitochondrialcontent and number of detected features

- "manual": Filters cells based on user-defined thresholds for
  nCount_RNA, nFeature_RNA, and percent_mt

- "mads": Filters cells based on median absolute deviations from the
  median for nCount_RNA, nFeature_RNA, and percent_mt
