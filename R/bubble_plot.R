bubblePlot = function(so,features,palette="RdBu",assay="SCT",ident="seurat_clusters"){
    Idents(so) = ident
    
    dotplot = DotPlot(so,features=features)
    pct.expression = select(dotplot$data,pct.exp,id,features.plot)
    
    avgExpression = AverageExpression(so,assay=assay)[[assay]]
    avgExpression = as.matrix(avgExpression)
    avgExpression.df = melt(avgExpression)
    colnames(avgExpression.df) = c("group","Gene","AvgExp")
    avgExpression.df$PctExp = 0
    for (i in 1:nrow(avgExpression.df)){
        avgExpression.df[i,4]=pct.expression$pct.exp[which(pct.expression$id==avgExpression.df[ident,i] & pct.expression$features.plot==avgExpression.df$Gene[i])]
    }

    plot = ggplot(avgExpression.df,aes(x=group,y=Gene,size=Pct.Exp,color=AvgExp))+
        geom_point()+
        scale_color_distiller(palette=palette)+
        theme_bw()+
        theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
    } 

    return (plot)

}