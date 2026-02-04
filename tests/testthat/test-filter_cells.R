test_that("filter_miQC", {
  pbmc <- suppressMessages(import_pbmc())
  set.seed(42)
  pbmc_filt <- filter_cells(pbmc,method = "miQC")
  expect_equal(length(which(pbmc_filt$keep == "keep")), 2140, tolerance = 5)
})

test_that("filter_manual", {
  pbmc <- suppressMessages(import_pbmc())
  set.seed(42)
  pbmc_filt <- filter_cells(pbmc,method = "manual")
  expect_equal(length(which(pbmc_filt$keep == "keep")), 2510, tolerance = 5)
})

test_that("filter_mads", {
  pbmc <- suppressMessages(import_pbmc())
  set.seed(42)
  pbmc_filt <- filter_cells(pbmc,method = "mads")
  expect_equal(length(which(pbmc_filt$keep == "keep")), 2590, tolerance = 5)
})
