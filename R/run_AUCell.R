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
  if (!inherits(so, "Seurat")) {
    rlang::abort("`so` must be a valid Seurat object.")
  }

  if (!is.list(gene_sets) || length(gene_sets) == 0) {
    rlang::abort("`gene_sets` must be a non-empty list of gene vectors.")
  }

  if (is.null(names(gene_sets)) || any(names(gene_sets) == "")) {
    rlang::abort("`gene_sets` must be a named list.")
  }

  if (!"SCT" %in% names(so@assays)) {
    rlang::abort("Seurat object must contain an `SCT` assay.")
  }

  expr_matrix <- SeuratObject::GetAssayData(
    so,
    assay = "SCT",
    layer = "data"
  )

  if (nrow(expr_matrix) == 0 || ncol(expr_matrix) == 0) {
    rlang::abort("Expression matrix is empty.")
  }

  gene_sets <- GSEABase::GeneSetCollection(
    lapply(
      seq_along(gene_sets),
      function(index) {
        GSEABase::GeneSet(
          geneIds = as.character(gene_sets[[index]]),
          setName = names(gene_sets)[index]
        )
      }
    )
  )

  # Find intersection of gene names between gene sets and Seurat object
  gene_sets <- AUCell::subsetGeneSets(gene_sets, rownames(expr_matrix))

  if (length(gene_sets) == 0) {
    rlang::abort(
      "No genes from `gene_sets` were found in the expression matrix."
    )
  }

  # Run AUCell rankings and calculations
  cells_rankings <- AUCell::AUCell_buildRankings(
    expr_matrix,
    nCores = 4,
    plotStats = FALSE
  )
  cells_AUC <- AUCell::AUCell_calcAUC(gene_sets, cells_rankings)

  cells_assignment <- AUCell::AUCell_exploreThresholds(
    cells_AUC,
    plotHist = FALSE,
    assignCells = TRUE
  )

  cells_assignment <- lapply(
    cells_assignment,
    function(single_assignment) {
      selected_cells <- single_assignment$assignment
      assignment_vector <- colnames(expr_matrix) %in% selected_cells
      names(assignment_vector) <- colnames(expr_matrix)
      single_assignment$assignment <- assignment_vector
      single_assignment
    }
  )
  # selectedThresholds <- getThresholdSelected(cells_assignment)

  # AUCell_plotTSNE(tSNE=obj@reductions$umap@cell.embeddings[,1:2],
  #  exprMat=exprMatrix,plots = c("AUC","binaryAUC"),#, "binaryAUC", "AUC","expression"),
  #  cellsAUC=cells_AUC[1:4,])

  aucell_output <- list(
    scores = as.matrix(AUCell::getAUC(cells_AUC)),
    assignment = cells_assignment
  )
  return(aucell_output)
}
