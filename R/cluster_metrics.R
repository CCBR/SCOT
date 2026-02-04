#' cluster_metrics: Calculates clustering metric scores
#'
#' @description
#' Uses the ClusterSim package to calculate metrics for cluster consistency.
#' These metrics do not require a priori knowledge of the truth set, which makes
#' them ideal for scRNA clustering.
#'
#' @details
#' Calculates three consistency scores to evaluate clustering effectiveness:
#' Calinski-Harabasz: Lower scores preferred.
#' Davies-Bouldin: Higher scores preferred.
#' Silhouette score: Higher scores preferred. Recommended to only run when data
#' has fewer than 90,000 data points
#'
#' @param so A Seurat single cell RNA (scRNA) object
#' @param cluster_list A vector of cluster resolution names (e.g.
#' "SCT_snn_res.0.2")
#' @param dims An integer vector of the principal components to be used in
#' calculating the cluster metrics. May not exceed the number of principal
#' components available
#' @param reduction A string for the reduction to be used in calculating scores
#' @param silhouette A Boolean for whether to include silhouette scoring
#'
#' @return Returns a table of clustering metrics for each cluster resolution
#' selected
#'
#' @export
#'
cluster_metrics <- function(
  so,
  cluster_list,
  dims = 1:20,
  reduction = "pca",
  silhouette = TRUE
) {
  embed_mat <- SeuratObject::Embeddings(so, reduction = reduction)[, dims]

  cluster_scores <- matrix(ncol = 2, nrow = length(cluster_list))
  colnames(cluster_scores)[1:2] <- c("CalinskiHarabasz", "DaviesBouldin")
  if (ncol(so) > 90000) {
    silhouette <- FALSE
  }
  if (silhouette == TRUE) {
    distance <- stats::dist(embed_mat)
    cluster_scores <- cbind(
      cluster_scores,
      vector(length = length(cluster_list))
    )
    colnames(cluster_scores)[3] <- "Silhouette"
  }
  rownames(cluster_scores) <- cluster_list

  for (i in seq_along(cluster_list)) {
    cluster_label <- cluster_list[i]
    clusters <- as.numeric(unlist(so[[cluster_label]]))

    ch_score <- clusterSim::index.G1(x = embed_mat, cl = clusters)
    db_score <- clusterSim::index.DB(x = embed_mat, cl = clusters)$DB
    cluster_scores[i, 1:2] <- c(ch_score, db_score)
    if (silhouette == TRUE) {
      sil_score <- clusterSim::index.S(d = distance, cl = clusters)
      cluster_scores[i, 3] <- sil_score
    }
  }

  return(cluster_scores)
}
