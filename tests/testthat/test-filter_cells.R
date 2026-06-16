test_that("filter_miQC", {
  set.seed(42)
  pbmc_filt <- filter_cells(pbmc_processed, method = "miQC")
  expect_equal(length(which(pbmc_filt$keep == "keep")), 2140, tolerance = 5)
})

test_that("filter_miQC_case", {
  set.seed(42)
  pbmc_filt <- filter_cells(pbmc_processed, method = "mIqC")
  expect_equal(length(which(pbmc_filt$keep == "keep")), 2140, tolerance = 5)
})

test_that("filter_manual", {
  set.seed(42)
  pbmc_filt <- filter_cells(pbmc_processed, method = "manual")
  expect_equal(length(which(pbmc_filt$keep == "keep")), 2510, tolerance = 5)
})

test_that("filter_manual_case", {
  set.seed(42)
  pbmc_filt <- filter_cells(pbmc_processed, method = "MAnUal")
  expect_equal(length(which(pbmc_filt$keep == "keep")), 2510, tolerance = 5)
})

test_that("filter_mads", {
  set.seed(42)
  pbmc_filt <- filter_cells(pbmc_processed, method = "mads")
  expect_equal(length(which(pbmc_filt$keep == "keep")), 2590, tolerance = 5)
})

test_that("filter_mads_case", {
  set.seed(42)
  pbmc_filt <- filter_cells(pbmc_processed, method = "maDS")
  expect_equal(length(which(pbmc_filt$keep == "keep")), 2590, tolerance = 5)
})

test_that("filter_nonexistent", {
  set.seed(42)
  expect_warning(
    pbmc_filt <- filter_cells(pbmc_processed, method = "whatever"),
    "No valid filtering"
  )
  expect_equal(ncol(pbmc_filt), 2638)
})
