#' splitFeaturePlot: Splits FeaturePlot by groups
#' 
#' @description
#' Adapts the FeaturePlot function from Seurat to split by groups and arrange by user specification
#' 
#' @param so A Seurat object
#' @param features A list of features (e.g. genes) to plot on a projection
#' @param splitIdent A metadata identity for splitting the image
#' @param label Boolean for whether to label the plot with the current Idents value of the Seurat object
#' @param ncol Integer value for number of columns in plot
#' @param nrow Integer value for number of rows in plot
#' @param min.cutoff Vector of minimum cutoff values for each feature, may specify quantile in the form of 'q##' where '##' is the quantile (eg, 1, 10)
#' @param max.cutoff Vector of maximum cutoff values for each feature, may specify quantile in the form of 'q##' where '##' is the quantile (eg, 1, 10)
#' @param plotImage Boolean for whether to export plots to current visualization window
#' @param returnList Boolean for whether to export list of plots
#' @param slot Data slot from which to extract values
#' @param order Boolean for whether to show highest feature values at the front of the image
#' @param reduction Dimensionality reduction to use
#' 
#' @return A list of ggplot2 plots


splitFeaturePlot = function(so,features,splitIdent,label=F,ncol=NA,nrow=NA,min.cutoff=NA,max.cutoff=NA,plotImage=F,returnList=F,slot="scale.data",order=F,reduction=NULL){

    plotList = list()
    plotOutput=list()
    if(is.null(reduction)){
        embed = so@reductions[[DefaultDimReduc(so)]]@cell.embeddings
    }else{
        embed = Embeddings(so, reduction=reduction)
    }

    for (feature in features){
        for (ident in unique(unlist(so[[splitIdent]]))){
            plotList[[feature]][[ident]] = FeaturePlot(so[,which(so[[splitIdent]] == ident)],features=feature,label=label,
                ,min.cutoff=min.cutoff,max.cutoff=max.cutoff,slot=slot,order=order,reduction=reduction) +
            xlim(range(embed[,1]))+ ylim (range(embed[,2]))+
            ggtitle(ident)
        }
    }
    for (feature in names(plotList)){
        if (is.na(ncol) & is.na(nrow)){ncol=length(plotList[[feature]]);nrow=1}
        if (is.na(ncol) & !is.na(nrow)){ncol = ceiling(length(plotList[[feature]])/nrow)}
        if (is.na(nrow) & !is.na(ncol)){nrow = ceiling(length(plotList[[feature]])/ncol)}

        plot.print = ggarrange(plotlist=plotList[[feature]], ncol=ncol,nrow=nrow,common.legend=T,legend="right")
        plot.print = annotate_figure(plot.print,top=text_grob(feature,face="bold",size=14))
        plotOutput[[feature]] = plot.print
        if(plotImage == T){
            dev.new()
            print(plot.print)        
        }
    }
    
    
    if(returnList==T){
        return(plotOutput)
    }

}
