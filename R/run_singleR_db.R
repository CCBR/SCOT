#' run_singleR_db: Runs SingleR with the built-in databases from celldex
#'
#' @description Uses the run_singleR function with the provided default
#' human or mouse references
#'
#' @param so_in A given Seurat single cell object
#' @param species Indicates human (hg19 or hg38) or mouse (mm10) references to
#' be used
#'
#'
#' @export
#'
#' @return A Seurat single cell object with predicted cell type annotations
run_singleR_db <- function(so_in, species) {
  # NOTE BEFORE NEXT PUSH: BUILD FUNCTION TO RETRIEVE SINGLE CELL REFERENCE OBJECTS EXTERNALLY

  # fetch_singleR_references <- function(ref = NULL){

  # <- celldex::HumanPrimaryCellAtlasData()
  # }

  # fetch_singleR_references() #Populates references into R environment

  cell_ont <- ontoProc::getOnto("cellOnto")
  if (species == "hg38" || species == "hg19") {
    so_in$HPCA_main <- run_singleR(
      so_in,
      fetch_celldex_ref("hpca"),
      "label.main"
    )
    so_in$HPCA_fine <- run_singleR(
      so_in,
      fetch_celldex_ref("hpca"),
      "label.fine"
    )
    so_in$HPCA_ont <- run_singleR(
      so_in,
      fetch_celldex_ref("hpca"),
      "label.ont"
    )
    so_in$HPCA_ont <- as.vector(cell_ont$name[so_in$HPCA_ont])

    so_in$BP_encode_main <- run_singleR(
      so_in,
      fetch_celldex_ref("BP_encode"),
      "label.main"
    )
    so_in$BP_encode_fine <- run_singleR(
      so_in,
      fetch_celldex_ref("BP_encode"),
      "label.fine"
    )
    so_in$BP_encode_ont <- run_singleR(
      so_in,
      fetch_celldex_ref("BP_encode"),
      "label.ont"
    )
    so_in$BP_encode_ont <- as.vector(cell_ont$name[so_in$BP_encode_ont])

    so_in$monaco_main <- run_singleR(
      so_in,
      fetch_celldex_ref("monaco"),
      "label.main"
    )
    so_in$monaco_fine <- run_singleR(
      so_in,
      fetch_celldex_ref("monaco"),
      "label.fine"
    )
    so_in$monaco_ont <- run_singleR(
      so_in,
      fetch_celldex_ref("monaco"),
      "label.ont"
    )
    so_in$monaco_ont <- as.vector(cell_ont$name[so_in$monaco_ont])

    so_in$immu_cell_exp_main <- run_singleR(
      so_in,
      fetch_celldex_ref("dice"),
      "label.main"
    )
    so_in$immu_cell_exp_fine <- run_singleR(
      so_in,
      fetch_celldex_ref("dice"),
      "label.fine"
    )
    so_in$immu_cell_exp_ont <- run_singleR(
      so_in,
      fetch_celldex_ref("dice"),
      "label.ont"
    )
    so_in$immu_cell_exp_ont <- as.vector(cell_ont$name[so_in$immu_cell_exp_ont])
    so_in$annot <- so_in$HPCA_main
  } else if (species == "mm10") {
    so_in$immgen_main <- run_singleR(
      so_in,
      fetch_celldex_ref("immgen"),
      "label.main"
    )
    so_in$immgen_fine <- run_singleR(
      so_in,
      fetch_celldex_ref("immgen"),
      "label.fine"
    )
    so_in$immgen_ont <- run_singleR(
      so_in,
      fetch_celldex_ref("immgen"),
      "label.ont"
    )
    so_in$immgen_ont <- as.vector(cell_ont$name[so_in$immgen_ont])

    so_in$mouseRNAseq_main <- run_singleR(
      so_in,
      fetch_celldex_ref("mouseRNAseq"),
      "label.main"
    )
    so_in$mouseRNAseq_fine <- run_singleR(
      so_in,
      fetch_celldex_ref("mouseRNAseq"),
      "label.fine"
    )
    so_in$mouseRNAseq_ont <- run_singleR(
      so_in,
      fetch_celldex_ref("mouseRNAseq"),
      "label.ont"
    )
    so_in$mouseRNAseq_ont <- as.vector(cell_ont$name[so_in$mouseRNAseq_ont])

    so_in$annot <- so_in$immgen_main
  }
  return(so_in)
}
