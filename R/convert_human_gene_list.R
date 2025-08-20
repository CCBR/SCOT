#' convert_human_gene_list
#'
#' @description Uses the ENSEMBL database to convert human gene names to mouse
#' gene names
#'
#' @param genes Character vector of human genes
#'
#' @import AnnotationDbi
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import Orthology.eg.db
#'
#' @export
#'
#' @return Character vector of mapped mouse genes


convert_human_gene_list <- function(genes) {
  egs <- mapIds(org.Hs.eg.db, gns, "ENTREZID", "SYMBOL")
  mapped <- select(Orthology.eg.db, egs, "Mus.musculus", "Homo.sapiens")
  mapped$MUS <- mapIds(
    org.Mm.eg.db,
    as.character(mapped$Mus.musculus),
    "SYMBOL", "ENTREZID"
  )
  return(as.character(unlist(mapped$MUS)))
}
