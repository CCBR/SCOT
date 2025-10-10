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
#'
#' @import Seurat
#'
#' @export
#'
#' @return A Seurat single cell object with normalized read counts

preprocess_sample <- function(so_in, species, npcs_in) {
  # assign genes depending on species input
  # TODO cc.genes is not initialized
  if (species == "hg38" || species == "hg19") {
    print("--proccesing human data")
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
  } else if (species == "mm10") {
    print("--proccesing mouse data")
    s.genes <- convert_human_gene_list(cc.genes$s.genes)
    g2m.genes <- convert_human_gene_list(cc.genes$g2m.genes)
  }

  # process
  so_1 <- NormalizeData(so_in,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    assay = "RNA"
  )
  so_2 <- ScaleData(so_1, assay = "RNA")
  so_3 <- CellCycleScoring(so_2,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE
  )
  so_4 <- SCTransform(so_3)
  so_out <- seurat_clustering(so_4, npcs_in)
  return(so_out)
}
