#' convert_mtx_to_h5: Converts the triple .mtx files from single cell
#' experiments to a .h5 counts file
#'
#' @description Uses existing tools to generate .h5 files, which can be easier
#' to import or transfer
#'
#' @details Older versions of CellRanger, in addition to some other single cell
#' alignment tools, such as PipSeq and DropSeq, produce a trio of files in an
#' output director that pointed to all the necessary features:
#' - matrix.mtx.gz: A matrix counts file with cells in the columns and features
#'  in the rows
#' - barcodes.tsv.gz: A tab-separated file with cell/droplet barcode sequences
#' - features.tsv.gz: A tab-separated file with gene names (or other features)
#'
#' @param sample_name A character string to be used to name the sample
#' @param mtx_file Path to the matrix.mtx.gz file
#' @param barcodes_file Path to the barcodes.tsv.gz file
#' @param features_file Path to the features.tsv.gz file
#'
#' @import DropletUtils
#' @import Seurat
#'
#' @export

convert_mtx_to_h5 <- function(
  sample_name,
  mtx_file,
  features_file,
  barcodes_file
) {
  counts <- Seurat::ReadMtx(
    mtx = mtx_file,
    cells = barcodes_file,
    features = features_file
  )
  outfile <- paste0(sample_name, ".h5")
  DropletUtils::write10xCounts(x = counts, path = outfile)
}
