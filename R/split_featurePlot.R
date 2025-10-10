#' split_featurePlot: Splits FeaturePlot by groups
#'
#' @description Adapts the FeaturePlot function from Seurat to split by groups
#' and arrange by user specification
#'
#' @param so A Seurat object
#' @param features A list of features (e.g. genes) to plot on a projection
#' @param split_ident A metadata identity for splitting the image
#' @param label Boolean for whether to label the plot with the current Idents
#' value of the Seurat object
#' @param ncol Integer value for number of columns in plot
#' @param nrow Integer value for number of rows in plot
#' @param min.cutoff Vector of minimum cutoff values for each feature, may
#' specify quantile in the form of 'q##' where '##' is the quantile (eg, 1, 10)
#' @param max.cutoff Vector of maximum cutoff values for each feature, may
#' specify quantile in the form of 'q##' where '##' is the quantile (eg, 1, 10)
#' @param plot_image Boolean for whether to export plots to current
#' visualization window
#' @param return_list Boolean for whether to export list of plots
#' @param slot Data slot from which to extract values
#' @param order Boolean for whether to show highest feature values at the front
#' of the image
#' @param reduction Dimensionality reduction to use
#'
#' @export
#'
#' @return A list of ggplot2 plots
split_featurePlot <- function(
    so,
    features,
    split_ident,
    label = FALSE,
    ncol = NA, nrow = NA,
    min.cutoff = NA, max.cutoff = NA,
    plot_image = FALSE,
    return_list = FALSE,
    slot = "scale.data",
    order = FALSE,
    reduction = NULL) {
  plot_list <- list()
  plot_output <- list()
  if (is.null(reduction)) {
    embed <- so@reductions[[SeuratObject::DefaultDimReduc(so)]]@cell.embeddings
  } else {
    embed <- Seurat::Embeddings(so, reduction = reduction)
  }

  # TODO refactor with map or lapply
  for (feature in features) {
    for (ident in unique(unlist(so[[split_ident]]))) {
      plot_list[[feature]][[ident]] <- Seurat::FeaturePlot(
        so[, which(so[[split_ident]] == ident)],
        features = feature, label = label,
        min.cutoff = min.cutoff, max.cutoff = max.cutoff,
        slot = slot, order = order, reduction = reduction
      ) +
        ggplot2::xlim(range(embed[, 1])) + ggplot2::ylim(range(embed[, 2])) +
        ggplot2::ggtitle(ident)
    }
  }
  # TODO refactor with map or lapply
  for (feature in names(plot_list)) {
    if (is.na(ncol) & is.na(nrow)) {
      ncol <- length(plot_list[[feature]])
      nrow <- 1
    }
    if (is.na(ncol) & !is.na(nrow)) {
      ncol <- ceiling(length(plot_list[[feature]]) / nrow)
    }
    if (is.na(nrow) & !is.na(ncol)) {
      nrow <- ceiling(length(plot_list[[feature]]) / ncol)
    }

    plot_print <- ggpubr::ggarrange(
      plotlist = plot_list[[feature]],
      ncol = ncol, nrow = nrow,
      common.legend = TRUE, legend = "right"
    )
    plot_print <- ggpubr::annotate_figure(plot_print,
      top = ggpubr::text_grob(feature, face = "bold", size = 14)
    )
    plot_output[[feature]] <- plot_print
    if (plot_image == TRUE) {
      grDevices::dev.new()
      print(plot_print)
    }
  }


  if (return_list == TRUE) {
    return(plot_output)
  }
}
