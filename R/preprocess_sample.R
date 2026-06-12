#' preprocess_sample: Performs an initial preprocess on a Seurat object
#'
#' @description Conducts preliminary normalization/scaling, cell cycle scoring,
#' and clustering on a Seurat single cell object
#'
#' @details The function performs the following processes on the data:assay
#' - Cell cycle scoring
#' - Normalization and scaling of RNA counts with log normalization
#' - Normalization and scaling of RNA counts with SCTransform
#' - Principal components analysis, neighbor identification, UMAP dimensionality
#' reduction and preliminary clustering (through `seurat_clustering` function)
#'
#' @param so_in A Seurat single cell object with RNA counts
#' @param species Species of the sample (hg38, hg19, or mm10)
#' @param npcs_in Number of principal components to use for dimensionality
#' reduction
#' @param replace Logical indicating whether to sample with replacement for
#' background gene selection in cell cycle scoring (default: FALSE). Note: this
#' parameter is retained for compatibility but currently has no effect, as
#' replacement is determined by the \code{ctrl} computation below.
#' @param num_ctrl Number of control genes to use for cell cycle scoring
#' background correction (default: 100)
#'
#'
#' @export
#'
#' @return A Seurat single cell object with normalized read counts
preprocess_sample <- function(
  so_in,
  species,
  npcs_in,
  replace = FALSE,
  num_ctrl = 100
) {
  if (is.null(so_in) || !inherits(so_in, "Seurat")) {
    stop("so_in must be a Seurat object")
  }
  if (!species %in% c("hg38", "hg19", "mm10")) {
    stop("species must be one of: hg38, hg19, mm10")
  }
  if (species == "hg38" || species == "hg19") {
    print("--proccesing human data")
    s.genes <- Seurat::cc.genes$s.genes
    g2m.genes <- Seurat::cc.genes$g2m.genes
  } else if (species == "mm10") {
    print("--proccesing mouse data")
    s.genes <- convert_human_gene_list(Seurat::cc.genes$s.genes)
    g2m.genes <- convert_human_gene_list(Seurat::cc.genes$g2m.genes)
  }

  # process
  so_1 <- Seurat::NormalizeData(
    so_in,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    assay = "RNA"
  )
  so_2 <- Seurat::ScaleData(so_1, assay = "RNA")
  #nbin <- 24L
  #ctrl <- min(round(nrow(so_2) * 0.25), nrow(so_2) %/% nbin - 1L)
  so_3 <- Seurat::CellCycleScoring(
    so_2,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE,
    ctrl = ,
    use.e = FALSE
  )
  so_4 <- Seurat::SCTransform(so_3)
  so_out <- seurat_clustering(so_4, npcs_in)
  return(so_out)
}
