# convert_mtx_to_h5: Converts the triple .mtx files from single cell experiments to a .h5 counts file

Uses existing tools to generate .h5 files, which can be easier to import
or transfer

## Usage

``` r
convert_mtx_to_h5(sample_name, mtx_file, features_file, barcodes_file)
```

## Arguments

- sample_name:

  A character string to be used to name the sample

- mtx_file:

  Path to the matrix.mtx.gz file

- features_file:

  Path to the features.tsv.gz file

- barcodes_file:

  Path to the barcodes.tsv.gz file

## Details

Older versions of CellRanger, in addition to some other single cell
alignment tools, such as PipSeq and DropSeq, produce a trio of files in
an output director that pointed to all the necessary features:

- matrix.mtx.gz: A matrix counts file with cells in the columns and
  features in the rows

- barcodes.tsv.gz: A tab-separated file with cell/droplet barcode
  sequences

- features.tsv.gz: A tab-separated file with gene names (or other
  features)
