# fetch_celldex_ref: Helper function to import annotation references from celldex

Retrieves one of the 7 available references in celldex. Should be
cache-agnostic by declaring the cache location as the local working
directory

## Usage

``` r
fetch_celldex_ref(ref_name, cache = ".")
```

## Arguments

- ref_name:

  Character string to be used for retrieving references from the celldex
  package

- cache:

  Character string to identify cache download location

## Value

A celldex single cell object in SCE data structure
