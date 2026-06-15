test_that("human singleR_db", {
  pbmc_annot <- run_singleR_db(pbmc_processed, species = "hg38")

  #Check annotations of HPCA subsets
  expect_equal(
    length(which(pbmc_annot$HPCA_main == "T_cells")),
    1394,
    tolerance = 5
  )
  expect_equal(
    length(grep("NK_cell", pbmc_annot$HPCA_fine)),
    156,
    tolerance = 5
  )
  expect_equal(
    length(which(pbmc_annot$HPCA_ont == "monocyte")),
    503,
    tolerance = 5
  )
  #Check annotations of BP ENCODE subsets
  expect_equal(
    length(which(pbmc_annot$BP_encode_main == "T-cells")),
    449,
    tolerance = 5
  )
  expect_equal(
    length(which(pbmc_annot$BP_encode_fine == "NK cells")),
    166,
    tolerance = 5
  )
  expect_equal(
    length(which(pbmc_annot$BP_encode_ont == "monocyte")),
    672,
    tolerance = 5
  )
  #Check annotations of Database of Immune Cell Expression subset
  expect_equal(
    length(which(pbmc_annot$immu_cell_exp_main == "T cells, CD4+")),
    929,
    tolerance = 5
  )
  expect_equal(
    length(which(pbmc_annot$immu_cell_exp_fine == "NK cells")),
    305,
    tolerance = 5
  )
  expect_equal(
    length(grep("monocyte", pbmc_annot$immu_cell_exp_ont)),
    671,
    tolerance = 5
  )
  #Check Monaco Immune annotation subsets
  expect_equal(
    length(which(pbmc_annot$monaco_main == "CD4+ T cells")),
    889,
    tolerance = 5
  )
  expect_equal(
    length(which(pbmc_annot$monaco_fine == "Natural killer cells")),
    147,
    tolerance = 5
  )
  expect_equal(
    length(grep("monocyte", pbmc_annot$monaco_ont)),
    619,
    tolerance = 5
  )
})

test_that("mouse singleR_db", {
  mouse_bm_annot <- run_singleR_db(mouse_bm_processed, species = "mm10")

  #Check ImmGen annotation subsets
  expect_equal(
    length(which(mouse_bm_annot$immgen_main == "Stem cells")),
    513,
    tolerance = 5
  )
  expect_equal(
    length(which(mouse_bm_annot$immgen_fine == "Stem cells (SC.MEP)")),
    288,
    tolerance = 5
  )
  expect_equal(
    length(which(
      mouse_bm_annot$immgen_ont == "megakaryocyte-erythroid progenitor cell"
    )),
    334,
    tolerance = 5
  )
  #Check MouseRNAseq annotations subsets
  expect_equal(
    length(which(mouse_bm_annot$mouseRNAseq_main == "Granulocytes")),
    454,
    tolerance = 5
  )
  expect_equal(
    length(which(mouse_bm_annot$mouseRNAseq_fine == "Erythrocytes)")),
    469,
    tolerance = 5
  )
  expect_equal(
    length(which(mouse_bm_annot$mouseRNAseq_ont == "granulocyte")),
    452,
    tolerance = 5
  )
})

test_that("Absent species", {
  expect_error(
    run_singleR_db(pbmc_processed, species = "mmul"),
    "No valid mouse or human genome or species label submitted"
  )
})

test_that("Missing input object", {
  expect_error(
    run_singleR_db(species = "mm10"),
    "argument .* is missing"
  )
})
