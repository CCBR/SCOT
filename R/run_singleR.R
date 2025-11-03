#' run_singleR: Runs SingleR directly on the Seurat object
#'
#' @description Wrapper script for the SingleR function
#'
#' @param so A Seurat single cell object
#' @param ref_file A reference file in a SingleCellExperiment format, such as those obtained from SingleR/CellDex package
#' @param label The label identity to be used. Must be a column header in the metadata of `ref_file`
#'
#' @export
#'
#' @return A vector of pruned cell type labels
run_singleR <- function(so, ref_file, label) {
  obj <- Seurat::DietSeurat(so, graphs = "umap")
  sce <- Seurat::as.SingleCellExperiment(so, assay = "SCT")
  s <- SingleR::SingleR(test = sce, ref = ref_file, labels = ref_file[[label]])
  return(s$pruned.labels)
}
