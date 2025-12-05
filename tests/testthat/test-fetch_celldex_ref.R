test_that("hpca_import", {
  hpca <- fetch_celldex_ref("hpca")
  hpca_2 <- fetch_celldex_ref("HumanPrimaryCellAtlasData")
  expect_equal(typeof(hpca), "S4")
  expect_equal(typeof(hpca_2), "S4")

  expect_equal(ncol(hpca), 713)
  expect_equal(ncol(hpca_2), 713)

  expect_equal(nrow(hpca), 19363)
  expect_equal(nrow(hpca_2), 19363)
})

test_that("bp_encode_import", {
  bp <- fetch_celldex_ref("blueprint_encode")
  bp_2 <- fetch_celldex_ref("BP_encode")
  bp_3 <- fetch_celldex_ref("BlueprintEncodeData")
  expect_equal(typeof(bp), "S4")
  expect_equal(typeof(bp_2), "S4")
  expect_equal(typeof(bp_3), "S4")

  expect_equal(ncol(bp), 259)
  expect_equal(ncol(bp_2), 259)
  expect_equal(ncol(bp_3), 259)

  expect_equal(nrow(bp), 19859)
  expect_equal(nrow(bp_2), 19859)
  expect_equal(nrow(bp_3), 19859)
})

test_that("monaco_import", {
  monaco <- fetch_celldex_ref("monaco")
  monaco_2 <- fetch_celldex_ref("MonacoImmuneData")
  expect_equal(typeof(monaco), "S4")
  expect_equal(typeof(monaco_2), "S4")

  expect_equal(ncol(monaco), 114)
  expect_equal(ncol(monaco_2), 114)

  expect_equal(nrow(monaco), 46077)
  expect_equal(nrow(monaco_2), 46077)
})

test_that("dice_import", {
  dice <- fetch_celldex_ref("dice")
  dice_2 <- fetch_celldex_ref("DatabaseImmuneCellExpressionData")
  expect_equal(typeof(dice), "S4")
  expect_equal(typeof(dice_2), "S4")

  expect_equal(ncol(dice), 1561)
  expect_equal(ncol(dice_2), 1561)

  expect_equal(nrow(dice), 48043)
  expect_equal(nrow(dice_2), 48043)
})
test_that("noversh_import", {
  noversh <- fetch_celldex_ref("hematopoietic")
  noversh_2 <- fetch_celldex_ref("NovershternHematopoieticData")
  noversh_3 <- fetch_celldex_ref("novershtern")
  expect_equal(typeof(noversh), "S4")
  expect_equal(typeof(noversh_2), "S4")
  expect_equal(typeof(noversh_3), "S4")

  expect_equal(ncol(noversh), 211)
  expect_equal(ncol(noversh_2), 211)
  expect_equal(ncol(noversh_3), 211)

  expect_equal(nrow(noversh), 13276)
  expect_equal(nrow(noversh_2), 13276)
  expect_equal(nrow(noversh_3), 13276)
})
test_that("immgen_import", {
  immgen <- fetch_celldex_ref("immgen")
  immgen_2 <- fetch_celldex_ref("ImmGenData")
  expect_equal(typeof(immgen), "S4")
  expect_equal(typeof(immgen_2), "S4")

  expect_equal(ncol(immgen), 830)
  expect_equal(ncol(immgen_2), 830)

  expect_equal(nrow(immgen), 22134)
  expect_equal(nrow(immgen_2), 22134)
})

test_that("mouseRNAseq_import", {
  mm_rna <- fetch_celldex_ref("mouseRNAseq")
  mm_rna_2 <- fetch_celldex_ref("MouseRNAseqData")
  expect_equal(typeof(mm_rna), "S4")
  expect_equal(typeof(mm_rna_2), "S4")

  expect_equal(ncol(mm_rna), 358)
  expect_equal(ncol(mm_rna_2), 358)

  expect_equal(nrow(mm_rna), 21214)
  expect_equal(nrow(mm_rna_2), 21214)
})

test_that("non_existent_ref", {
  ref_test <- fetch_celldex_ref("blech")
  expect_null(ref_test)
  expect_error(fetch_celldex_ref(1))
  expect_error(fetch_celldex_ref())
})
