#' add_VDJ_metadata: Appends V(D)J alignments from CellRanger to an existing
#' Seurat object
#' 
#' @description Wrapper for djvdj to assign V(D)J alignments to cells in a 
#' Seurat object. Can be used for TCRa/TCRb, TCRg/TCRd, and BCR
#' 
#' @param so A Seurat single cell RNA object
#' @param vdj_file A path to a V(D)J annotation from CellRanger
#' 
#' @import djvdj
#' @import Seurat
#' 
#' @export 
#' 
#'@return A Seurat object with V(D)J metadata

add_VDJ_metadata = function(so, vdj_file){
  prefix <- basename(vdj_file)

  #Used in the case of multiple V(D)J alignments
  if (prefix == "outs"){
    prefix <- NA
  } else {
    prefix <- paste0(prefix, "_")
  }

  so <- djvdj::import_vdj(
    input = so,
    vdj_dir = vdj_file,
    filter_paired = FALSE,
    define_clonotypes = "cdr3_gene",
    prefix = prefix
  )

  return(so)
}