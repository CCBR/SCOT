test_that("object_size", {
  pbmc <- suppressMessages(import_pbmc())
  expect_equal(ncol(pbmc), 2638)
})

test_that("cell_type_annotation", {
  pbmc <- suppressMessages(import_pbmc())
  set.seed(42)
  ref <- fetch_celldex_ref("hpca")

  singleR_output <- run_singleR(pbmc, ref, "label.main")

  expect_equal(length(which(singleR_output == "T_cells")), 1394)
  expect_equal(length(which(singleR_output == "B_cell")), 335)
  expect_equal(length(which(singleR_output == "Monocyte")), 619)
  expect_equal(length(which(singleR_output == "NK_cell")), 185)
})
