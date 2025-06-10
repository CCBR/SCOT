#' clusterMetrics: Calculates clustering metric scores
#' 
#' @description 
#' Uses the ClusterSim package to calculate metrics for cluster consistency. These metrics do not require a priori knowledge of the truth set, which makes them ideal for scRNA clustering.
#' 
#' @details 
#' Calculates three consistency scores to evaluate clustering effectiveness:
#' Calinski-Harabasz: Lower scores preferred.
#' Davies-Bouldin: Higher scores preferred.
#' Silhouette score: Higher scores preferred. Recommended to only run when data has fewer than 90,000 data points
#' 
#' @param so A Seurat single cell RNA (scRNA) object
#' @param clusterList A vector of cluster resolution names (e.g. "SCT_snn_res.0.2") 
#' @param dims An integer vector of the principal components to be used in calculating the cluster metrics. May not exceed the number of princiapl components available
#' @param reduction A string for the reduction to be used in calcualting scores
#' @param silhouette A Boolean for whether to include silhouette scoring
#' 
#' @return Returns a table of clustering metrics for each cluster resolution selected

clusterMetrics = function(so, clusterList, dims=1:20, reduction="pca",silhouette=T){
    embedMat = Embeddings(so,reduction=reduction)[,dims]

    clusterScores = matrix(ncol=2,nrow=length(clusterList))
    colnames(clusterScores)[1:2] = c("CalinskiHarabasz","DaviesBouldin")
    if (ncol(so)>90000){silhouette=F}
    if(silhouette==T){
        distance = dist(embedMat)
        clusterScores = cbind(clusterScores,vector(length=length(clusterList)))
        colnames(clusterScores)[3]="Silhouette"
    }
    rownames(clusterScores) = clusterList

    for (i in 1:length(clusterList)){
        clusterLabel = clusterList[i]
        clusters=as.numeric(unlist(so[[clusterLabel]]))

        ch_score = clusterSim::index.G1(x=embedMat,cl=clusters)
        db_score = clusterSim::index.DB(x=embedMat,cl=clusters)$DB
        clusterScores[i,1:2] = c(ch_score,db_score)
        if(silhouette==T){
            sil_score=clusterSim::index.S(d=distance,cl=clusters)
            clusterScores[i,3]=sil_score
        }

    }


    return(clusterScores)
}