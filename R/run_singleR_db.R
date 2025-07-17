#' run_singleR_db: Runs SingleR with the built-in databases from celldex
#'
#' @description Uses the run_singleR function with the provided default human or mouse references
#'
#' @param so_in A given Seurat single cell object
#' @param species Indicates human (hg19 or hg38) or mouse (mm10) references to be used
#'
#' @return A Seurat single cell object with predicted cell type annotations

run_singleR_db <- function(so_in, species) {
  cell_ont <- ontoProc::getOnto("cellOnto")
  if (species == "hg38" || species == "hg19") {
    so_in$HPCA_main <- run_singleR(
      so_in, celldex::HumanPrimaryCellAtlasData(), "label.main"
    )
    so_in$HPCA_fine <- run_singleR(
      so_in, celldex::HumanPrimaryCellAtlasData(), "label.fine"
    )
    so_in$HPCA_ont <- run_singleR(
      so_in, celldex::HumanPrimaryCellAtlasData(), "label.ont"
    )
    so_in$HPCA_ont <- cellOnt$name[so_in$HPCA_ont]

    so_in$BP_encode_main <- run_singleR(
      so_in, celldex::BlueprintEncodeData(), "label.main"
    )
    so_in$BP_encode_fine <- run_singleR(
      so_in, celldex::BlueprintEncodeData(), "label.fine"
    )
    so_in$BP_encode_ont <- run_singleR(
      so_in, celldex::BlueprintEncodeData()(), "label.ont"
    )
    so_in$BP_encode_ont <- cellOnt$name[so_in$BP_encode_ont]

    so_in$monaco_main <- run_singleR(
      so_in, celldex::MonacoImmuneData(), "label.main"
    )
    so_in$monaco_fine <- run_singleR(
      so_in, celldex::MonacoImmuneData(), "label.fine"
    )
    so_in$monaco_ont <- run_singleR(
      so_in, celldex::MonacoImmuneData(), "label.ont"
    )
    so_in$monaco_ont <- cellOnt$name[so_in$monaco_ont]

    so_in$immu_cell_exp_main <- run_singleR(
      so_in, celldex::DatabaseImmuneCellExpressionData(),
      "label.main"
    )
    so_in$immu_cell_exp_fine <- run_singleR(
      so_in, celldex::DatabaseImmuneCellExpressionData(),
      "label.fine"
    )
    so_in$immu_cell_exp_ont <- run_singleR(
      so_in, celldex::DatabaseImmuneCellExpressionData(), "label.ont"
    )
    so_in$immu_cell_exp_ont <- cellOnt$name[so_in$immu_cell_exp_ont]
    so_in$annot <- so_in$HPCA_main
  } else if (species == "mm10") {
    so_in$immgen_main <- run_singleR(so_in, celldex::ImmGenData(), "label.main")
    so_in$immgen_fine <- run_singleR(so_in, celldex::ImmGenData(), "label.fine")
    so_in$immgen_ont <- run_singleR(so_in, celldex::ImmGenData(), "label.ont")
    so_in$immgen_ont <- cell_ont$names[so_in$immgen_ont]

    so_in$mouseRNAseq_main <- run_singleR(
      so_in, celldex::MouseRNAseqData(), "label.main"
    )
    so_in$mouseRNAseq_fine <- run_singleR(
      so_in, celldex::MouseRNAseqData(), "label.fine"
    )
    so_in$mouseRNAseq_ont <- run_singleR(
      so_in, celldex::MouseRNAseqData(), "label.ont"
    )
    so_in$mouseRNAseq_ont <- cellOnt$name[so_in$mouseRNAseq_ont]

    so_in$annot <- so_in$immgen_main
  }
  return(so_in)
}
