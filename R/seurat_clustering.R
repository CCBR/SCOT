#' seurat_clustering: A quick end-to-end clustering of a Seurat object
#'
#' @description Wrapper to perform PCA, neighbor-finding, UMAP calculations and clustering
#'
#' @param so_in A Seurat single cell object
#' @param npcs_in Number of principal components to be used to calculate
#' neighbors, UMAP projection, and clustering
#' @param resolution A resolution value for clustering. Higher values will
#' produce more clusters (Default 0.8)
#' @param algorithm A value of 1 (Louvain), 2 (Leiden), or 3 (Slow Local Moving
#' (SLM)) to select the clustering algorithm (Default 3)
#'
#' @export
#'
#' @return The updated Seurat single cell object with recalculated PCA,
#' neighbors, UMAP and clusters
seurat_clustering <- function(so_in, npcs_in, resolution = 0.8, algorithm = 3) {
  so <- Seurat::RunPCA(
    object = so_in,
    features = Seurat::VariableFeatures(object = so_in),
    verbose = FALSE,
    npcs = 50
  )
  so <- Seurat::FindNeighbors(so, dims = 1:npcs_in)
  so <- Seurat::FindClusters(
    so,
    resolution = 0.8,
    algorithm = 3,
    verbose = TRUE
  )
  so <- Seurat::RunUMAP(so, dims = 1:npcs_in, n.components = 3)
  return(so)
}
