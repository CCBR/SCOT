run_AUCell <- function(so, gene_sets) {
  set.seed(42)

  # Retrieve normalized counts matrix
  expr <- FetchData(
    so,
    vars = rownames(so@assays$SCT@data), slot = "SCT", layer = "data"
  )
  expr_matrix <- t(expr)
  # Find intersection of gene names between gene sets and Seurat object
  gene_sets <- subsetGeneSets(gene_sets, rownames(expr_matrix))
  gene_sets <- setGeneSetNames(gene_sets, newNames = paste(names(gene_sets), " (", nGenes(gene_sets), "g)", sep = ""))

  # Run AUCell rankings and calculations
  cells_rankings <- AUCell_buildRankings(expr_matrix, nCores = 4, plotStats = FALSE)
  cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings)

  cells_assignment <- AUCell_exploreThresholds(
    cells_AUC,
    plotHist = FALSE, assign = TRUE
  )
  # selectedThresholds <- getThresholdSelected(cells_assignment)

  # AUCell_plotTSNE(tSNE=obj@reductions$umap@cell.embeddings[,1:2],
  #  exprMat=exprMatrix,plots = c("AUC","binaryAUC"),#, "binaryAUC", "AUC","expression"),
  #  cellsAUC=cells_AUC[1:4,])

  aucell_output <- list()

  aucell_output$scores <- cells_AUC
  aucell_output$assignment <- cells_assignment

  return(aucell_output)
}
