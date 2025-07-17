#' filter_doublets: Identifies and removes doublets or multiplets in scRNA
#' samples
#'
#' @description
#' Wrapper for DoubletFinder and scDblFinder for identifying doublets
#'
#' @details
#' Runs the latest version of DoubletFinder and/or scDblFinder and returns
#' a filtered scRNA object with doublets removed
#'
#' @param so_in A provided Seurat scRNA object
#' @param doublet_finder_method Character string indicating the use of
#' "DoubletFinder" (default), "scDblFinder", "consensus" removal of doublets,
#' or the "union" of total doublets identified in both algorithms
#'
#' @return

filter_doublets <- function(so_in, doublet_finder_method = "DoubletFinder") {
  if (doublet_finder_method == "DoubletFinder") {
    sweep_res_list <- paramSweep(so_in, PCs = 1:10, sct = T)
    sweep_stats <- summarizeSweep(sweep_res_list, GT = FALSE)
    bcmvn <- find.pK(sweep_stats)

    ## Homotypic Doublet Proportion Estimate
    homotypic_prop <- modelHomotypic(so_in$annot)
    perc <- 0.005 * (length(colnames(so_in)) / 1000)
    nExp_poi <- round(perc * length(colnames(so_in))) # dfso@cell.names
    nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))

    ## Run DoubletFinder with varying classification stringencies
    dfso <- doubletFinder(so_in,
      pN = 0.25, pK = 0.09,
      nExp = nExp_poi,
      reuse.pANN = NULL, PCs = 1:10, sct = TRUE
    )

    pAAN <- tail(names(dfso@meta.data), 2)[1]
    dfso <- doubletFinder(dfso,
      pN = 0.25, pK = 0.09,
      nExp = nExp_poi_adj,
      reuse.pANN = pAAN, PCs = 1:10, sct = T
    )
    so_in$doubletFinder_label <- dfso[[tail(names(dfso@meta.data), 1)]]
    so_in <- subset(
      so_in,
      cells = colnames(so_in)[which(so_in$doubletFinder_label == "Singlet")]
    )
  } else if (doublet_finder_method == "scDblFinder") {
    sce <- as.SingleCellExperiment(so_in)
    set.seed(42)
    sce_dbl <- scDblFinder(sce) %>% suppressWarnings()
    so_in$scDblFinder_label <- sce_dbl$scDblFinder.class
    so_in <- subset(
      so_in,
      cells = colnames(so_in)[which(so_in$scDblFinder_label == "singlet")]
    )
  } else if (doublet_finder_method == "union") {
    sweep_res_list <- paramSweep(so_in, PCs = 1:10, sct = T)
    sweep_stats <- summarizeSweep(sweep_res_list, GT = FALSE)
    bcmvn <- find.pK(sweep_stats)

    ## Homotypic Doublet Proportion Estimate
    homotypic_prop <- modelHomotypic(so_in$annot)
    perc <- 0.005 * (length(colnames(so_in)) / 1000)
    nExp_poi <- round(perc * length(colnames(so_in))) # dfso@cell.names
    nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))

    ## Run DoubletFinder with varying classification stringencies
    dfso <- doubletFinder(so_in,
      pN = 0.25, pK = 0.09,
      nExp = nExp_poi,
      reuse.pANN = NULL, PCs = 1:10, sct = TRUE
    )

    pAAN <- tail(names(dfso@meta.data), 2)[1]
    dfso <- doubletFinder(dfso,
      pN = 0.25, pK = 0.09,
      nExp = nExp_poi_adj,
      reuse.pANN = pAAN, PCs = 1:10, sct = T
    )
    so_in$doubletFinder_label <- dfso[[tail(names(dfso@meta.data), 1)]]

    sce <- as.SingleCellExperiment(so_in)
    set.seed(42)
    sce_dbl <- scDblFinder(sce) %>% suppressWarnings()
    so_in$scDblFinder_label <- sce_dbl$scDblFinder.class
    so_in <- subset(
      so_in,
      cells = colnames(so_in)[intersect(
        which(so_in$doubletFinder_label == "Singlet"),
        which(so_in$scDblFinder_label == "singlet")
      )]
    )
  } else if (doublet_finder_method == "consensus") {
    sweep_res_list <- paramSweep(so_in, PCs = 1:10, sct = T)
    sweep_stats <- summarizeSweep(sweep_res_list, GT = FALSE)
    bcmvn <- find.pK(sweep_stats)

    ## Homotypic Doublet Proportion Estimate
    homotypic_prop <- modelHomotypic(so_in$annot)
    perc <- 0.005 * (length(colnames(so_in)) / 1000)
    nExp_poi <- round(perc * length(colnames(so_in))) # dfso@cell.names
    nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))

    ## Run DoubletFinder with varying classification stringencies
    dfso <- doubletFinder(so_in,
      pN = 0.25, pK = 0.09,
      nExp = nExp_poi,
      reuse.pANN = NULL, PCs = 1:10, sct = TRUE
    )

    pAAN <- tail(names(dfso@meta.data), 2)[1]
    dfso <- doubletFinder(dfso,
      pN = 0.25, pK = 0.09,
      nExp = nExp_poi_adj,
      reuse.pANN = pAAN, PCs = 1:10, sct = T
    )
    so_in$doubletFinder_label <- dfso[[tail(names(dfso@meta.data), 1)]]

    sce <- as.SingleCellExperiment(so_in)
    set.seed(42)
    sce_dbl <- scDblFinder(sce) %>% suppressWarnings()
    so_in$scDblFinder_label <- sce_dbl$scDblFinder.class
    so_in <- subset(
      so_in,
      cells = colnames(so_in)[unique(c(
        which(so_in$doubletFinder_label == "Singlet"),
        which(so_in$scDblFinder_label == "singlet")
      ))]
    )
  } else {
    cat("No valid method selected. Returning initial object")
  }

  return(so_in)
}
