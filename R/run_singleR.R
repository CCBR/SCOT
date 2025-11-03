#' run_singleR: Runs SingleR directly on the Seurat object
#'
#' @description Wrapper script for the SingleR function
#'
#' @param so A Seurat single cell object
#' @param ref_file A reference file in a SingleCellExperiment format, such as those obtained from SingleR/CellDex package
#' @param label The label identity to be used. Must be a column header in the metadata of `ref_file`
#'
#' @import SingleR
#' @import Seurat
#'
#' @export
#'
#' @return A vector of pruned cell type labels


run_singleR <- function(so, ref_file, label) {
  obj <- DietSeurat(so, graphs = "umap")
  sce <- as.SingleCellExperiment(so, assay = "SCT")
  ref <- ref_file
  s <- SingleR(test = sce, ref = ref, labels = ref[[label]])
  return(s$pruned.labels)
}
