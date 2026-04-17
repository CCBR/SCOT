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
#' @param de_table A differential expression table generated from the Seurat
#' FindMarkers function. Contains the column headers `p_val`, `avg_log2FC`,
#' `pct.1`, `pct.2`, and `p_val_adj`
#'
#' @export
#'
#' @return Returns a named vector of GSEA scores for genes
calc_sc_gsea_score <- function(de_table) {
  gsea_score_vect <- sign(de_table$avg_log2FC) *
    -log10(de_table$p_val) *
    pmax(de_table$pct.1, de_table$pct.2) +
    de_table$avg_log2FC

  # pmax compares pct.1 and pct.2 for each gene and returns the maximum value
  gsea_score_vect[which(de_table$p_val == 0)] <- (sign(de_table$avg_log2FC) *
    500 *
    pmax(de_table$pct.1, de_table$pct.2) +
    de_table$avg_log2FC)[which(de_table$p_val == 0)]

  names(gsea_score_vect) <- rownames(de_table)

  return(gsea_score_vect)
}
