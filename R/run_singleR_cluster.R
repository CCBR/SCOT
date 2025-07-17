#' run_singleR_cluster: Runs SingleR on average expression of clusters
#'
#' @description Wrapper function to run SingleR cell type annotation on the
#' current default identities of a Seurat scRNA object, e.g. clusters
#'
#' @param so_in A Seurat single cell object
#' @param ref_file A SingleR compatible single cell reference object
#' (SingleCellExperiment object)
#' @param label A cell type label from the metadata column headers
#'
#' @return A character vector of matched cell type annotations based on
#' clusters
#'

run_singleR_cluster <- function(so_in, ref_file, label) {
  avg <- AverageExpression(so_in, assays = "SCT")
  avg <- as.data.frame(avg)
  ref <- refFile
  s <- SingleR(test = as.matrix(avg), ref = ref, labels = ref[[label]])

  clust_annot <- s$labels
  names(clust_annot) <- colnames(avg)
  names(clust_annot) <- gsub("SCT.", "", names(clust_annot))

  annot_vect <- clust_annot[match(so_in$seurat_clusters, names(clust_annot))]
  names(annot_vect) <- colnames(so_in)
  return(annot_vect)
}
