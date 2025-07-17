runAUCell = function(so,geneSets){
    set.seed(42)
  #expr = FetchData(obj,vars=rownames(obj@assays$RNA@data) ,slot="data")
  #Retrieve normalized counts matrix
  expr = FetchData(so,vars=rownames(so@assays$SCT@data),slot="SCT", layer="data")
  exprMatrix <- t(expr)
  #Find intersection of gene names between gene sets and Seurat object
  geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix))
  geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets), " (", nGenes(geneSets) ,"g)", sep=""))
  countsPerGene <- apply(exprMatrix, 1, function(x) sum(x>0))
  #Run AUCell rankings and calculations
  cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=4, plotStats=TRUE)
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
  cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, assign=TRUE)
  selectedThresholds <- getThresholdSelected(cells_assignment)
  #AUCell_plotTSNE(tSNE=obj@reductions$umap@cell.embeddings[,1:2],
  #  exprMat=exprMatrix,plots = c("AUC","binaryAUC"),#, "binaryAUC", "AUC","expression"),
  #  cellsAUC=cells_AUC[1:4,])
  return(cells_AUC)    
}