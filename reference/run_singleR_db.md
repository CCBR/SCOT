# run_singleR_db: Runs SingleR with the built-in databases from celldex

Uses the run_singleR function with the provided default human or mouse
references

## Usage

``` r
run_singleR_db(so_in, species)
```

## Arguments

- so_in:

  A given Seurat single cell object

- species:

  Indicates human (hg19 or hg38) or mouse (mm10) references to be used

## Value

A Seurat single cell object with predicted cell type annotations
