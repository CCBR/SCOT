#' run_AUCell: A wrapper script for running AUCell for per cell gene set
#' enrichment analysis
#'
#' @description Runs AUCell on basic inputs
#'
#' @details AUCell evaluates the likelihood of a gene set being enriched in an
#' individual cells using the Area Under the Curve (AUC). The wrapper script
#' uses a Seurat object and a list of gene sets to calculate the AUC for each
#' gene set and also provides the assignment likelihood, based on the
#' distribution of the AUCs across the dataset.
#'
#' @param so A Seurat object
#' @param gene_sets A list of gene sets, such as those obtained via GSEABase
#'
#' @export
#'
#' @return A list object containing the AUCell scores for each gene set and the
#' assignment of the cells as to the likelihood of the cells overexpressing the
#' gene set. These values can be added to the original Seurat object ex post
#' facto.

run_AUCell <- function(so, gene_sets) {
  set.seed(42)

  # Retrieve normalized counts matrix
  expr <- Seurat::FetchData(
    so,
    vars = rownames(so@assays$SCT@data), slot = "SCT", layer = "data"
  )
  expr_matrix <- t(expr)
  # Find intersection of gene names between gene sets and Seurat object
  gene_sets <- AUCell::subsetGeneSets(gene_sets, rownames(expr_matrix))
  gene_sets <- AUCell::setGeneSetNames(gene_sets, newNames = paste(names(gene_sets), " (", AUCell::nGenes(gene_sets), "g)", sep = ""))

  # Run AUCell rankings and calculations
  cells_rankings <- AUCell::AUCell_buildRankings(expr_matrix, nCores = 4, plotStats = FALSE)
  cells_AUC <- AUCell::AUCell_calcAUC(gene_sets, cells_rankings)

  cells_assignment <- AUCell::AUCell_exploreThresholds(
    cells_AUC,
    plotHist = FALSE, assignCells = TRUE
  )
  # selectedThresholds <- getThresholdSelected(cells_assignment)

  # AUCell_plotTSNE(tSNE=obj@reductions$umap@cell.embeddings[,1:2],
  #  exprMat=exprMatrix,plots = c("AUC","binaryAUC"),#, "binaryAUC", "AUC","expression"),
  #  cellsAUC=cells_AUC[1:4,])

  aucell_output <- list(
    scores = cells_AUC,
    assignment = cells_assignment
  )
  return(aucell_output)
}
