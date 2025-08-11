#' bubble_plot: Flexible dotplot
#' 
#' @description 
#' 
#' @param 
#' 
#' @return 
#' 


bubble_plot <- function(
  so, features, palette = "RdBu", assay = "SCT", ident = "seurat_clusters"
) {
  Idents(so) <- ident

  dotplot <- DotPlot(so, features = features)
  pct_expression <- select(dotplot$data, pct.exp, id, features.plot)

  avg_expression <- AverageExpression(so, assay = assay)[[assay]]
  avg_expression <- as.matrix(avg_expression)
  avg_expression_df <- melt(avg_expression)
  colnames(avg_expression_df) <- c("group", "Gene", "AvgExp")
  avg_expression_df$PctExp <- 0
  for (i in seq_len(avg_expression_df)) {
    avg_expression_df[i, 4] <- pct_expression$pct.exp[which(
      pct_expression$id == avg_expression_df[ident, i] &
        pct_expression$features.plot == avg_expression_df$Gene[i]
    )]
  }

  plot <- ggplot(avg_expression_df, aes(
    x = group, y = Gene, size = Pct.Exp, color = AvgExp
  )) +
    geom_point() +
    scale_color_distiller(palette = palette) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

  return(plot)
}
