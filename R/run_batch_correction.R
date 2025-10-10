#' run_batch_correction: Performs batch correction within Seurat v5 framework
#'
#' @description
#' Performs batch correction using Seurat v5, with additional dimensionality
#' reduction (UMAP), clustering, and cluster-based cell type annotation.
#'
#' @details
#' Utilizes the Seurat v5 framework to streamline batch correction using one of
#' the following options:
#' - SCVI
#' - LIGER
#' - Canonical Correlation Analysis (CCA, or Seurat "integrate")
#' - Reciprocal PCA (RPCA)
#' - Harmony
#' Follows up with re-running dimensionality reduction with PCA, neighbor
#' finding, and UMAP projection. Also runs clustering with the slow local
#' moving algorithm at a series of resolutions selected by the user. Final
#' step conducts cell type annotation with SingleR using the average gene
#' expression vector of each cluster
#'
#' @param so_in A merged Seurat object containing all samples prior to batch
#' correction
#' @param npcs The number of principal components to use for dimensionality
#' reduction and neighbor identification
#' @param species "hg19", "hg38", or "mm10". Used to determine databases used
#' for cell type annotation
#' @param resolution_list A vector of resolutions ranging from 0 to 2.0 for
#' clustering. Smaller resolutions produce fewer clusters (Default)
#' @param method_in A character string to indicate which batch correction method
#' to use
#' @param reduction_in A character string for naming the PCA produced
#' @param vars_to_regress A character vector for variables to regress out when
#' running SCTransform normalization and batch correction
#' @param conda_env A character string for indicating which conda environment
#' contains the necessary packages
#'
#' @import Seurat
#' @import SeuratWrappers
#' @import harmony
#' @import rliger
#'
#' @export
#'
#' @return A batch corrected Seurat object with updated UMAP projections and
#' clusters
#'
run_batch_correction <- function(
    so_in, npcs, species, resolution_list, method_in,
    reduction_in = c(0.2, 0.4, 0.5, 0.8, 1.0),
    vars_to_regress = NULL, conda_env = "") {
  # data variables must be initialized to silence the R CMD check note:
  #    'no visible binding for global variable'
  v_list <- NULL

  # TODO recommend using package::function syntax instead of importing entire packages

  # set assay to RNA to avoid double transform/normalization
  DefaultAssay(so_in) <- "RNA"

  # integration method for
  ### SCVI
  ### LIGER
  ### harmony,rpca,cca
  if (method_in == "scVIIntegration") {
    message("--running SCVI integration")

    so_transform <- NormalizeData(so_in)
    so_variable <- FindVariableFeatures(so_transform)
    so_scaled <- ScaleData(so_variable)
    so_pca <- RunPCA(so_scaled)

    so_integrate <- IntegrateLayers(
      object = so_pca, method = scVIIntegration,
      new.reduction = "integrated.scvi",
      conda_env = conda_env, dims = 1:npcs
    )
  } else if (method_in == "LIGER") {
    message("--running LIGER")

    # preprocess
    so_norm <- NormalizeData(so_in)
    so_norm <- FindVariableFeatures(so_norm)
    so_norm <- ScaleData(so_norm, do.center = FALSE)
    so_norm <- RunOptimizeALS(so_norm, k = npcs, lambda = 5)
    so_integrate <- RunQuantileNorm(so_norm)
  } else {
    message("--running SCT")

    # vars.to.regress is NULL by default
    so_transform <- SCTransform(so_in, vars.to.regress = v_list)

    # runPCA
    so_pca <- RunPCA(so_transform)

    so_integrate <- IntegrateLayers(
      object = so_pca, method = get(method_in),
      normalization.method = "SCT",
      verbose = FALSE, new.reduction = reduction_in
    )
  }

  # run neighbors, clusters
  so <- FindNeighbors(so_integrate, reduction = reduction_in, dims = 1:npcs)
  for (res in resolution_list) {
    so <- FindClusters(so, resolution = res, algorithm = 3)
  }

  # reduction
  so <- RunUMAP(so, reduction = reduction_in, dims = 1:npcs)

  # add cluster-based annotations
  cell_ont <- ontoProc::getOnto("cellOnto")
  if (species == "hg38" || species == "hg19") {
    so$clustAnnot_HPCA_main <- run_singleR_cluster(
      so, celldex::HumanPrimaryCellAtlasData(), "label.main"
    )
    so$clustAnnot_HPCA_fine <- run_singleR_cluster(
      so, celldex::HumanPrimaryCellAtlasData(), "label.fine"
    )
    so$clustAnnot_HPCA_ont <- run_singleR_cluster(
      so, celldex::HumanPrimaryCellAtlasData(), "label.ont"
    )
    so$clustAnnot_HPCA_ont <- cell_ont$name[so$clustAnnot_HPCA_ont]

    so$clustAnnot_BP_encode_main <- run_singleR_cluster(
      so, celldex::BlueprintEncodeData(), "label.main"
    )
    so$clustAnnot_BP_encode_fine <- run_singleR_cluster(
      so, celldex::BlueprintEncodeData(), "label.fine"
    )
    so$clustAnnot_BP_encode_ont <- run_singleR_cluster(
      so, celldex::BlueprintEncodeData(), "label.ont"
    )
    so$clustAnnot_BP_encode_ont <- cell_ont$name[so$clustAnnot_BP_encode_ont]
    so$clustAnnot_monaco_main <- run_singleR_cluster(
      so, celldex::MonacoImmuneData(), "label.main"
    )
    so$clustAnnot_monaco_fine <- run_singleR_cluster(
      so, celldex::MonacoImmuneData(), "label.fine"
    )
    so$clustAnnot_monaco_ont <- run_singleR_cluster(
      so, celldex::MonacoImmuneData(), "label.ont"
    )
    so$clustAnnot_monaco_ont <- cell_ont$name[so$clustAnnot_monaco_ont]
    so$clustAnnot_immu_cell_exp_main <- run_singleR_cluster(
      so, celldex::DatabaseImmuneCellExpressionData(), "label.main"
    )
    so$clustAnnot_immu_cell_exp_fine <- run_singleR_cluster(
      so, celldex::DatabaseImmuneCellExpressionData(), "label.fine"
    )
    so$clustAnnot_immu_cell_exp_ont <- run_singleR_cluster(
      so, celldex::DatabaseImmuneCellExpressionData(), "label.ont"
    )
    so$clustAnnot_immu_cell_exp_ont <- cell_ont$name[so$clustAnnot_immu_cell_exp_ont]
  } else if (species == "mm10") {
    so$clustAnnot_immgen_main <- run_singleR_cluster(
      so, celldex::ImmGenData(), "label.main"
    )
    so$clustAnnot_immgen_fine <- run_singleR_cluster(
      so, celldex::ImmGenData(), "label.fine"
    )
    so$clustAnnot_immgen_ont <- run_singleR_cluster(
      so, celldex::ImmGenData(), "label.ont"
    )
    so$clustAnnot_immgen_ont <- cell_ont$name[so$clustAnnot_immgen_ont]
    so$clustAnnot_mouseRNAseq_main <- run_singleR_cluster(
      so, celldex::MouseRNAseqData(), "label.main"
    )
    so$clustAnnot_mouseRNAseq_fine <- run_singleR_cluster(
      so, celldex::MouseRNAseqData(), "label.fine"
    )
    so$clustAnnot_mouseRNAseq_ont <- run_singleR_cluster(
      so, celldex::MouseRNAseqData(), "label.ont"
    )
    so$clustAnnot_mouseRNAseq_ont <- cell_ont$name[so$clustAnnot_mouseRNAseq_ont]
  }
  return(so)
}
