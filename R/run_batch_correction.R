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
#' @param reduction_in A character string for naming the reduction produced by
#' batch correction (default: "integrated")
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
  resolution_list = c(0.5, 1.0),
  method_in,
  reduction_in = "integrated",
  vars_to_regress = NULL,
  conda_env = ""
) {
  # data variables must be initialized to silence the R CMD check note:
  #    'no visible binding for global variable'
  v_list <- NULL

  if (is.null(so_in) || !inherits(so_in, "Seurat")) {
    stop("so_in must be a non-NULL Seurat object")
  }

  valid_species <- c("hg38", "hg19", "mm10")
  if (!species %in% valid_species) {
    stop(paste(
      "species not supported:",
      species,
      "- must be one of:",
      paste(valid_species, collapse = ", ")
    ))
  }

  valid_methods <- c(
    "scVIIntegration",
    "LIGER",
    "RPCAIntegration",
    "CCAIntegration",
    "HarmonyIntegration"
  )
  if (!method_in %in% valid_methods) {
    stop(paste(
      "method_in not supported:",
      method_in,
      "- must be one of:",
      paste(valid_methods, collapse = ", ")
    ))
  }

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
    # split RNA layers by batch before SCTransform (Seurat v5 integration workflow)
    split_col <- if ("batch" %in% colnames(so_in@meta.data)) {
      so_in$batch
    } else {
      SeuratObject::Idents(so_in)
    }
    so_split <- so_in
    so_split[["RNA"]] <- split(so_split[["RNA"]], f = split_col)
    so_transform <- Seurat::SCTransform(so_split, vars.to.regress = v_list)

    # runPCA
    so_pca <- Seurat::RunPCA(so_transform)

    # resolve method string to function from Seurat namespace
    method_fn <- tryCatch(
      utils::getFromNamespace(method_in, "Seurat"),
      error = function(e) stop(paste("method_in not supported:", method_in))
    )

    # k.weight must be < smallest batch size
    batch_sizes <- table(split_col)
    k_weight <- min(30L, min(batch_sizes) - 1L)
    so_integrate <- Seurat::IntegrateLayers(
      object = so_pca,
      method = method_fn,
      normalization.method = "SCT",
      verbose = FALSE,
      new.reduction = reduction_in,
      k.weight = k_weight
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
    .safe_annotate <- function(so, ref_name, label) {
      tryCatch(
        run_singleR_cluster(so, fetch_celldex_ref(ref_name), label),
        error = function(e) {
          warning(
            "run_singleR_cluster failed for ref=",
            ref_name,
            ", label=",
            label,
            ": ",
            conditionMessage(e)
          )
          rep(NA_character_, ncol(so))
        }
      )
    }

    so$clustAnnot_HPCA_main <- .safe_annotate(so, "hpca", "label.main")
    so$clustAnnot_HPCA_fine <- .safe_annotate(so, "hpca", "label.fine")
    so$clustAnnot_HPCA_ont <- .safe_annotate(so, "hpca", "label.ont")
    ont_vals <- so$clustAnnot_HPCA_ont
    so$clustAnnot_HPCA_ont <- ifelse(
      is.na(ont_vals),
      NA_character_,
      cell_ont$name[ont_vals]
    )

    so$clustAnnot_BP_encode_main <- .safe_annotate(
      so,
      "BP_encode",
      "label.main"
    )
    so$clustAnnot_BP_encode_fine <- .safe_annotate(
      so,
      "BP_encode",
      "label.fine"
    )
    so$clustAnnot_BP_encode_ont <- .safe_annotate(so, "BP_encode", "label.ont")
    ont_vals <- so$clustAnnot_BP_encode_ont
    so$clustAnnot_BP_encode_ont <- ifelse(
      is.na(ont_vals),
      NA_character_,
      cell_ont$name[ont_vals]
    )

    so$clustAnnot_monaco_main <- .safe_annotate(so, "monaco", "label.main")
    so$clustAnnot_monaco_fine <- .safe_annotate(so, "monaco", "label.fine")
    so$clustAnnot_monaco_ont <- .safe_annotate(so, "monaco", "label.ont")
    ont_vals <- so$clustAnnot_monaco_ont
    so$clustAnnot_monaco_ont <- ifelse(
      is.na(ont_vals),
      NA_character_,
      cell_ont$name[ont_vals]
    )

    so$clustAnnot_immu_cell_exp_main <- .safe_annotate(so, "dice", "label.main")
    so$clustAnnot_immu_cell_exp_fine <- .safe_annotate(so, "dice", "label.fine")
    so$clustAnnot_immu_cell_exp_ont <- .safe_annotate(so, "dice", "label.ont")
    ont_vals <- so$clustAnnot_immu_cell_exp_ont
    so$clustAnnot_immu_cell_exp_ont <- ifelse(
      is.na(ont_vals),
      NA_character_,
      cell_ont$name[ont_vals]
    )
  } else if (species == "mm10") {
    .safe_annotate <- function(so, ref_name, label) {
      tryCatch(
        run_singleR_cluster(so, fetch_celldex_ref(ref_name), label),
        error = function(e) {
          warning(
            "run_singleR_cluster failed for ref=",
            ref_name,
            ", label=",
            label,
            ": ",
            conditionMessage(e)
          )
          rep(NA_character_, ncol(so))
        }
      )
    }

    so$clustAnnot_immgen_main <- .safe_annotate(so, "immgen", "label.main")
    so$clustAnnot_immgen_fine <- .safe_annotate(so, "immgen", "label.fine")
    so$clustAnnot_immgen_ont <- .safe_annotate(so, "immgen", "label.ont")
    ont_vals <- so$clustAnnot_immgen_ont
    so$clustAnnot_immgen_ont <- ifelse(
      is.na(ont_vals),
      NA_character_,
      cell_ont$name[ont_vals]
    )

    so$clustAnnot_mouseRNAseq_main <- .safe_annotate(
      so,
      "mouseRNAseq",
      "label.main"
    )
    so$clustAnnot_mouseRNAseq_fine <- .safe_annotate(
      so,
      "mouseRNAseq",
      "label.fine"
    )
    so$clustAnnot_mouseRNAseq_ont <- .safe_annotate(
      so,
      "mouseRNAseq",
      "label.ont"
    )
    ont_vals <- so$clustAnnot_mouseRNAseq_ont
    so$clustAnnot_mouseRNAseq_ont <- ifelse(
      is.na(ont_vals),
      NA_character_,
      cell_ont$name[ont_vals]
    )
  }
  return(so)
}
