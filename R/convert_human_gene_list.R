#' convert_human_gene_list
#'
#' @description Uses the ENSEMBL database to convert human gene names to mouse
#' gene names
#'
#' @param genes Character vector of human genes
#'
#' @export
#'
#' @return Character vector of mapped mouse genes
convert_human_gene_list <- function(genes) {
  # data variables must be initialized to silence the R CMD check note:
  #    'no visible binding for global variable'
  gns <- NULL
  gns <- genes

  # TODO: make this function generic enough to convert any genome to any other
  egs <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, gns, "ENTREZID", "SYMBOL")
  mapped <- AnnotationDbi::select(
    Orthology.eg.db::Orthology.eg.db,
    egs,
    "Mus.musculus",
    "Homo.sapiens"
  )
  mapped$MUS <- AnnotationDbi::mapIds(
    org.Mm.eg.db::org.Mm.eg.db,
    as.character(mapped$Mus.musculus),
    "SYMBOL",
    "ENTREZID"
  )
  return(as.character(unlist(mapped$MUS)))
}
