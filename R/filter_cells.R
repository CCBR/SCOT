#' filter_cells: A wrapper for multiple scRNA filtering methods
#'
#' @description Wrapper function to filter cells from a Seurat object using
#' method of choice. Returns a Seurat object with filtering results
#'
#' @param so Seurat object to filter
#' @param method Filtering method to use: "miQC" (default), "manual", or "mads"
#' @param mads Number of median absolute deviations to use for "mads" method (default 3)
#' @param nCount_RNA_max Maximum nCount_RNA threshold for "manual" method (default 500000)
#' @param nCount_RNA_min Minimum nCount_RNA threshold for "manual" method (default 1000)
#' @param nFeature_RNA_max Maximum nFeature_RNA threshold for "manual" method (default 5000)
#' @param nFeature_RNA_min Minimum nFeature_RNA threshold for "manual" method (default 200)
#' @param percent_mt_max Maximum percent_mt threshold for "manual" method (default 10)
#' @param percent_mt_min Minimum percent_mt threshold for "manual" method (default 0)
#'
#' @details This function filters cells from a Seurat object using one of three
#' methods:
#' - "miQC": Uses the miQC package to filter cells based on a probabilistic
#' model based on mitochondrialcontent and number of detected features
#' - "manual": Filters cells based on user-defined thresholds for nCount_RNA,
#' nFeature_RNA, and percent_mt
#' - "mads": Filters cells based on median absolute deviations from the median for
#' nCount_RNA, nFeature_RNA, and percent_mt
#'
#' @return Returns a Seurat object with an additional metadata column "keep" indicating
#' whether each cell is to be kept or discarded
#'
#' @export
#'
filter_cells <- function(
  so,
  method = "miQC",
  mads = 3,
  nCount_RNA_max = 50000,
  nCount_RNA_min = 1000,
  nFeature_RNA_max = 5000,
  nFeature_RNA_min = 200,
  percent_mt_max = 10,
  percent_mt_min = 0
) {
  nFeature_RNA <- NULL
  so_filt <- BiocGenerics::subset(so, subset = nFeature_RNA > 200) # standard low cell read filter
  if (method == "miQC") {
    print("--filtering cells with miQC")

    so_qc <- SeuratWrappers::RunMiQC(
      so_filt,
      percent.mt = "percent.mt",
      nFeature_RNA = "nFeature_RNA",
      posterior.cutoff = 0.7,
      model.slot = "flexmix_model"
    )
    so_qc$keep <- so_qc$miQC.keep
  } else if (method == "manual") {
    so_qc <- so_filt
    so_qc$keep <- "discard"
    so_qc$keep[which(
      so_qc$nCount_RNA <= nCount_RNA_max &
        so_qc$nCount_RNA >= nCount_RNA_min &
        so_qc$nFeature_RNA <= nFeature_RNA_max &
        so_qc$nFeature_RNA >= nFeature_RNA_min &
        so_qc$percent.mt <= percent_mt_max &
        so_qc$percent.mt >= percent_mt_min
    )] <- "keep"
  } else if (method == "mads") {
    so_qc <- so_filt
    nCount_out <- Routliers::outliers_mad(
      so_qc$nCount_RNA,
      threshold = mads
    )$LL_CI_MAD
    nFeature_out <- Routliers::outliers_mad(
      so_qc$nFeature_RNA,
      threshold = mads
    )$LL_CI_MAD
    percent_mt_out <- Routliers::outliers_mad(
      so_qc$percent.mt,
      threshold = mads
    )$UL_CI_MAD
    so_qc$keep <- "discard"
    so_qc$keep[which(
      so_qc$nCount_RNA >= nCount_out &
        so_qc$nFeature_RNA >= nFeature_out &
        so_qc$percent.mt <= percent_mt_out
    )] <- "keep"
  }

  return(so_qc)
}
