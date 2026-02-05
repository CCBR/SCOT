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
#' @export
#'
#' @return A batch corrected Seurat object with updated UMAP projections and
#' clusters
#'
run_batch_correction <- function(
  so_in,
  npcs,
  species,
  resolution_list,
  method_in,
  reduction_in = c(0.2, 0.4, 0.5, 0.8, 1.0),
  vars_to_regress = NULL,
  conda_env = ""
) {
  # data variables must be initialized to silence the R CMD check note:
  #    'no visible binding for global variable'
  v_list <- NULL

  # TODO recommend using package::function syntax instead of importing entire packages

  # set assay to RNA to avoid double transform/normalization
  SeuratObject::DefaultAssay(so_in) <- "RNA"

  # integration method for
  ### SCVI
  ### LIGER
  ### harmony,rpca,cca
  if (method_in == "scVIIntegration") {
    message("--running SCVI integration")

    so_transform <- Seurat::NormalizeData(so_in)
    so_variable <- Seurat::FindVariableFeatures(so_transform)
    so_scaled <- Seurat::ScaleData(so_variable)
    so_pca <- Seurat::RunPCA(so_scaled)

    so_integrate <- Seurat::IntegrateLayers(
      object = so_pca,
      method = SeuratWrappers::scVIIntegration,
      new.reduction = "integrated.scvi",
      conda_env = conda_env,
      dims = 1:npcs
    )
  } else if (method_in == "LIGER") {
    message("--running LIGER")

    # preprocess
    so_norm <- Seurat::NormalizeData(so_in)
    so_norm <- Seurat::FindVariableFeatures(so_norm)
    so_norm <- Seurat::ScaleData(so_norm, do.center = FALSE)
    so_norm <- SeuratWrappers::RunOptimizeALS(so_norm, k = npcs, lambda = 5)
    so_integrate <- SeuratWrappers::RunQuantileNorm(so_norm)
  } else {
    message("--running SCT")

    # vars.to.regress is NULL by default
    so_transform <- Seurat::SCTransform(so_in, vars.to.regress = v_list)

    # runPCA
    so_pca <- Seurat::RunPCA(so_transform)

    so_integrate <- Seurat::IntegrateLayers(
      object = so_pca,
      method = get(method_in),
      normalization.method = "SCT",
      verbose = FALSE,
      new.reduction = reduction_in
    )
  }

  # run neighbors, clusters
  so <- Seurat::FindNeighbors(
    so_integrate,
    reduction = reduction_in,
    dims = 1:npcs
  )
  for (res in resolution_list) {
    so <- Seurat::FindClusters(so, resolution = res, algorithm = 3)
  }

  # reduction
  so <- Seurat::RunUMAP(so, reduction = reduction_in, dims = 1:npcs)

  # add cluster-based annotations
  cell_ont <- ontoProc::getOnto("cellOnto")
  if (species == "hg38" || species == "hg19") {
    so$clustAnnot_HPCA_main <- run_singleR_cluster(
      so,
      fetch_celldex_ref("hpca"),
      "label.main"
    )
    so$clustAnnot_HPCA_fine <- run_singleR_cluster(
      so,
      fetch_celldex_ref("hpca"),
      "label.fine"
    )
    so$clustAnnot_HPCA_ont <- run_singleR_cluster(
      so,
      fetch_celldex_ref("hpca"),
      "label.ont"
    )
    so$clustAnnot_HPCA_ont <- cell_ont$name[so$clustAnnot_HPCA_ont]

    so$clustAnnot_BP_encode_main <- run_singleR_cluster(
      so,
      fetch_celldex_ref("BP_encode"),
      "label.main"
    )
    so$clustAnnot_BP_encode_fine <- run_singleR_cluster(
      so,
      fetch_celldex_ref("BP_encode"),
      "label.fine"
    )
    so$clustAnnot_BP_encode_ont <- run_singleR_cluster(
      so,
      fetch_celldex_ref("BP_encode"),
      "label.ont"
    )
    so$clustAnnot_BP_encode_ont <- cell_ont$name[so$clustAnnot_BP_encode_ont]
    so$clustAnnot_monaco_main <- run_singleR_cluster(
      so,
      fetch_celldex_ref("monaco"),
      "label.main"
    )
    so$clustAnnot_monaco_fine <- run_singleR_cluster(
      so,
      fetch_celldex_ref("monaco"),
      "label.fine"
    )
    so$clustAnnot_monaco_ont <- run_singleR_cluster(
      so,
      fetch_celldex_ref("monaco"),
      "label.ont"
    )
    so$clustAnnot_monaco_ont <- cell_ont$name[so$clustAnnot_monaco_ont]
    so$clustAnnot_immu_cell_exp_main <- run_singleR_cluster(
      so,
      fetch_celldex_ref("dice"),
      "label.main"
    )
    so$clustAnnot_immu_cell_exp_fine <- run_singleR_cluster(
      so,
      fetch_celldex_ref("dice"),
      "label.fine"
    )
    so$clustAnnot_immu_cell_exp_ont <- run_singleR_cluster(
      so,
      fetch_celldex_ref("dice"),
      "label.ont"
    )
    so$clustAnnot_immu_cell_exp_ont <- cell_ont$name[
      so$clustAnnot_immu_cell_exp_ont
    ]
  } else if (species == "mm10") {
    so$clustAnnot_immgen_main <- run_singleR_cluster(
      so,
      fetch_celldex_ref("immgen"),
      "label.main"
    )
    so$clustAnnot_immgen_fine <- run_singleR_cluster(
      so,
      fetch_celldex_ref("immgen"),
      "label.fine"
    )
    so$clustAnnot_immgen_ont <- run_singleR_cluster(
      so,
      fetch_celldex_ref("immgen"),
      "label.ont"
    )
    so$clustAnnot_immgen_ont <- cell_ont$name[so$clustAnnot_immgen_ont]
    so$clustAnnot_mouseRNAseq_main <- run_singleR_cluster(
      so,
      fetch_celldex_ref("mouseRNAseq"),
      "label.main"
    )
    so$clustAnnot_mouseRNAseq_fine <- run_singleR_cluster(
      so,
      fetch_celldex_ref("mouseRNAseq"),
      "label.fine"
    )
    so$clustAnnot_mouseRNAseq_ont <- run_singleR_cluster(
      so,
      fetch_celldex_ref("mouseRNAseq"),
      "label.ont"
    )
    so$clustAnnot_mouseRNAseq_ont <- cell_ont$name[
      so$clustAnnot_mouseRNAseq_ont
    ]
  }
  return(so)
}
