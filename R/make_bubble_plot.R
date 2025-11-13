#' make a flexible bubble plot with Seurat::DotPlot
#'
#' @description Creates a bubble plot that shows average expression percent
#' expression for genes in single cell categories
#'
#' @param so A Seurat single cell RNA object
#' @param features A character vector of genes to plot
#' @param palette A color palette from ggplot2
#' @param assay The counts assay to use for determining expression
#' @param ident The categorical identity to classify groups of cells
#'
#' @export
#'
#' @return A ggplot2 figure
#'
make_bubble_plot <- function(so,
                             features,
                             palette = "RdBu",
                             assay = "SCT",
                             ident = "seurat_clusters") {
  # data variables must be initialized to silence the R CMD check note:
  #    'no visible binding for global variable'
  pct.exp <- id <- features.plot <- group <- Gene <- Pct.Exp <- AvgExp <- NULL

  Seurat::Idents(so) <- ident

  dotplot <- Seurat::DotPlot(so, features = features)
  pct_expression <- dplyr::select(dotplot$data, pct.exp, id, features.plot)

  avg_expression <- Seurat::AverageExpression(so, assays = assay)[[assay]][features,]
  avg_expression <- as.matrix(avg_expression)
  avg_expression_df <- reshape2::melt(avg_expression)
  
  # group and Gene appears to be in the wrong order
  ## original
  ## colnames(avg_expression_df) <- c("group", "Gene", "AvgExp")
  
  ## new
  colnames(avg_expression_df) <- c("Gene","group","AvgExp")
  avg_expression_df$PctExp <- 0

  # TODO refactor this with map
  for (i in 1:nrow(avg_expression_df)) {
    avg_expression_df[i, 4] <- pct_expression$pct.exp[which(
      pct_expression$id == avg_expression_df$group[i] &
        pct_expression$features.plot == avg_expression_df$Gene[i]
    )]
  }

  bubble_plot <- ggplot2::ggplot(
    avg_expression_df,
    ggplot2::aes(
      x = group,
      y = Gene,
      size = Pct.Exp,
      color = AvgExp
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_distiller(palette = palette) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ))

  return(bubble_plot)
}
