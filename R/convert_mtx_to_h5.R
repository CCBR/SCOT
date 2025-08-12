#' convert_mtx_to_h5: Converts the triple .mtx files from single cell experiments to a .h5 counts file
#'
#' @description
#'
#' @details
#'
#' @param
#'
#' @return

library(Seurat)
library(DropletUtils)

convert_mtx_to_h5 <- function(
    sample_name,
    mtx_file,
    features_file,
    barcodes_file) {
  counts <- Seurat::ReadMtx(mtx = mtx_file, cells = barcodes_file, features = features_file)
  outfile <- paste0(sample_name, ".h5")
  DropletUtils::write10xCounts(x = counts, path = outfile)
}
