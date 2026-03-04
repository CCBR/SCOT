test_that("statistics_added", {
  set.seed(42)
  pbmc_test <- calc_prelim_stats(load_plain_pbmc())
  expect_equal(mean(pbmc_test$percent.mt), 2.2, tolerance = 0.1)
  expect_equal(mean(pbmc_test$percent.RP), 37.3, tolerance = 0.2)
})
