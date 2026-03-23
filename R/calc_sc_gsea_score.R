#' calc_sc_gsea_score: Calculates scores for pre-ranked GSEA
#'
#' @description Creates a vector of scores that can be used to rank genes in
#' order to be submitted to GSEA
#'
#' @details The scoring metric uses a Seurat FindMarkers output table and
#' incorporates the p-value and the fold change to determine magnitude, while
#' using the percent expression values as weighting factors.  The basic formula
#' is:
#'
#' -log_10(p-value) * sign (log2FC) * max(pct.1, pct.2) + log2FC
#'
#' If the p-value is recorded as zero (i.e. below machine error), the formula is
#' adjusted
#'
#' 500 * sign(log2FC) * max(pct.1, pct.2) + log2FC
#'
#' This formula should preserve the ranking as indicated by magnitude of fold
#' change,while also diminishing the significance of lowly-expressed genes.
#'
#' @param deTable A differential expression table generated from the Seurat
#' FindMarkers function. Contains the column headers `p_val`, `avg_log2FC`,
#' `pct.1`, `pct.2`, and `p_val_adj`
#'
#' @export
#'
#' @return Returns a named vector of GSEA scores for genes
calc_sc_gsea_score <- function(deTable) {
  gseaScoreVect <- sign(deTable$avg_log2FC) *
    -log10(deTable$p_val) *
    pmax(deTable$pct.1, deTable$pct.2) +
    deTable$avg_log2FC

  gseaScoreVect[which(deTable$p_val == 0)] <- (sign(deTable$avg_log2FC) *
    500 *
    pmax(deTable$pct.1, deTable$pct.2) +
    deTable$avg_log2FC)[which(deTable$p_val == 0)]

  names(gseaScoreVect) <- rownames(deTable)

  return(gseaScoreVect)
}
