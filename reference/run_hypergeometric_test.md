# run_hypergeometric_test: Calculates the p-value of an overrepresentation test using the hypergeometric test

Calculates the likelihood that a selected reference gene set is
overrepresented in a provided gene list. This is a wrapper for the
`phyper` function in R that uses vectors of genes rather than requiring
users to understand the parameters of `phyper`. The probability is
evaluated as P(X \> x),where X is a potential value, and x is the actual
number of intersected genes.

## Usage

``` r
run_hypergeometric_test(selectedVect, refVect, worldSize, lowerTail = FALSE)
```

## Arguments

- selectedVect:

  A vector of selected genes from a dataset

- refVect:

  A vector of genes for a specific gene set to be evaluated for
  significance

- worldSize:

  An integer for the number of genes that are present in the dataset as
  a whole

- lowerTail:

  Calculates the p-value as likelihood of being in the upper tail
  (default) or the lower tail of the distribution

## Value

Returns a vector with the relevant statistics used to generate the
hypergeometric test
