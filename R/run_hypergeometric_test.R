#' run_hypergeometric_test: Calculates the p-value of an overrepresentation test
#' using the hypergeometric test
#'
#' @description Calculates the likelihood that a selected reference gene set is
#' overrepresented in a provided gene list. This is a wrapper for the `phyper`
#' function in R that uses vectors of genes rather than requiring users to
#' understand the parameters of `phyper`. The probability is evaluated as
#' P(X > x),where X is a potential value, and x is the actual number of
#' intersected genes.
#'
#' @param selectedVect A vector of selected genes from a dataset
#' @param refVect A vector of genes for a specific gene set to be evaluated for
#' significance
#' @param worldSize An integer for the number of genes that are present in the
#' dataset as a whole
#' @param lowerTail Calculates the p-value as likelihood of being in the upper
#' tail (default) or the lower tail of the distribution
#'
#' @import stats
#'
#' @export
#'
#' @return Returns a vector with the relevant statistics used to generate the
#' hypergeometric test
run_hypergeometric_test <- function(selectedVect,
                                    refVect,
                                    worldSize,
                                    lowerTail = FALSE) {
  intersectList <- intersect(selectedVect, refVect)
  inter.size <- length(intersectList)
  ref.size <- length(refVect)
  nonref.size <- worldSize - ref.size
  selection.size <- length(selectedVect)

  pVal <- 1
  if (inter.size > 0) {
    pVal <- phyper(
      q = inter.size,
      m = ref.size,
      n = nonref.size,
      k = selection.size,
      lower.tail = lowerTail
    )
  }

  statVector <- c(
    pVal,
    inter.size,
    selection.size,
    ref.size,
    worldSize,
    paste(intersectList, collapse = ";")
  )

  names(statVector) <- c(
    "P-value",
    "IntersectSize",
    "SelectGenesSize",
    "ReferenceSize",
    "WorldSize",
    "IntersectedGenes"
  )

  return(statVector)
}
