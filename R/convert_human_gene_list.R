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
convert_human_gene_list <- function(genes = NULL) {
  result <- character(0)

  # Validate input: must be character vector
  if (!is.character(genes)) {
    stop(paste0(
      "Input 'genes' must be a character vector, not ",
      paste(class(genes), collapse = " ")
    ))
  } else if (length(genes) == 0L) {
    warning("Input 'genes' is empty; returning character(0)")
  } else {
    # Map human SYMBOL to ENTREZID
    egs <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      genes,
      "ENTREZID",
      "SYMBOL"
    )

    # Map human ENTREZID to mouse ENTREZID
    mapped <- AnnotationDbi::select(
      Orthology.eg.db::Orthology.eg.db,
      as.character(egs),
      "Mus.musculus",
      "Homo.sapiens"
    )

    # Map mouse ENTREZID to SYMBOL
    mapped$MUS <- AnnotationDbi::mapIds(
      org.Mm.eg.db::org.Mm.eg.db,
      as.character(mapped$Mus.musculus),
      "SYMBOL",
      "ENTREZID"
    )

    result <- as.character(unlist(mapped$MUS))
  }

  return(result)
}
