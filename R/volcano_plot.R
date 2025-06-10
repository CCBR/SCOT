volcano_plot=function(deTable){
    log10_p = -log10(deTable$p_val_adj)
    log10_p[which(deTable$p_val_adj==0)] = 500
    avg_log2FC = deTable$avg_log2FC
    significance = vector(length=length(log10_p))
    significance[] = "NotSignificant"
    significance[which(abs(deTable$avg_log2FC) > 1.5)] = "avg_log2FC > 1.5"
    significance[which(deTable$p_val_adj<0.05)] = "p_val_adj < 0.05"
    significance[which(abs(deTable$avg_log2FC) > 1.5 & deTable$p_val_adj<0.05)] = "p_val_adj < 0.05 & avg_log2FC > 1.5"

    df = as.data.frame(cbind(avg_log2FC, log10_p, significance))
    rownames(df) = rownames(deTable)
    df[,1] = as.numeric(df[,1])
    df[,2] = as.numeric(df[,2])
    df[,3] = factor(df[,3],levels=c("NotSignificant","avg_log2FC > 1.5", "p_val_adj < 0.05", "p_val_adj < 0.05 & avg_log2FC > 1.5"))

    volcano = ggplot(df,aes(x=avg_log2FC,y=log10_p,col=significance))+
        geom_point()

}