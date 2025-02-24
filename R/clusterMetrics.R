#' clusterMetrics: Calculates clustering metric scores
#' 
#' @description 
#' Uses the ClusterSim package to calculate metrics for cluster consistency. These metrics do not 
#' require a priori knowledge of the truth set, which makes them ideal for scRNA clustering.
#' 
#' @details 
#' 
#' 
#' @param so 
clusterMetrics = function(so, clusterList, dims=1:20, reduction="pca",silhouette=T){
    embedMat = Embeddings(so,reduction=reduction)[,dims]

    clusterScores = matrix(ncol=2,nrow=length(clusterList))
    colnames(clusterScores)[1:2] = c("CalinskiHarabasz","DaviesBouldin")
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