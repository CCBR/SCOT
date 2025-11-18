#' fetch_celldex_ref: Helper function to import annotation references from
#' celldex
#'
#' @description Retrieves one of the 7 available references in celldex. Should
#' be cache-agnostic by declaring the cache location as the local working
#' directory
#'
#' @param ref_name Character string to be used for retrieving references from
#' the celldex package
#'
#' @export
#'
#' @return A celldex single cell object in SCE data structure
fetch_celldex_ref <- function(ref_name, cache = '.') {
  ref <- switch(ref_name,
    "hpca" = ,
    "HumanPrimaryCellAtlasData" = celldex::fetchReference(
      "hpca",
      version = "2024-02-26",
      realize.assays = TRUE, cache = "./"
    ),
    "blueprint_encode" = ,
    "BP_encode" = ,
    "bpencode" = ,
    "bp_encode" = ,
    "BlueprintEncodeData" = celldex::fetchReference(
      "blueprint_encode", "2024-02-26",
      realize.assays = TRUE, cache = "./"
    ),
    "monaco" = ,
    "MonacoImmuneData" = celldex::fetchReference(
      "monaco_immune", "2024-02-26",
      realize.assays = TRUE, cache = "./"
    ),
    "immu_cell_exp" = ,
    "DatabaseImmuneCellExpressionData" = ,
    "dice" = celldex::fetchReference(
      "dice", "2024-02-26",
      realize.assays = TRUE, cache = "./"
    ),
    "hematopoietic" = ,
    "NovershternHematopoieticData" = ,
    "novershtern_hematopoietic" = ,
    "novershtern" = celldex::fetchReference(
      "novershtern_hematopoietic", "2024-02-26",
      realize.assays = TRUE, cache = "./"
    ),
    "immgen" = ,
    "ImmGenData" = celldex::fetchReference(
      "immgen", "2024-02-26",
      realize.assays = TRUE, cache = "./"
    ),
    "mouseRNAseq" = ,
    "MouseRNAseqData" = celldex::fetchReference(
      "mouse_rnaseq", "2024-02-26",
      realize.assays = TRUE, cache = "./"
    )
  )
  return(ref)
}
