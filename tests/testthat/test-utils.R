test_that("abort_packages_not_installed works", {
  expect_error(abort_packages_not_installed("NOT_A_PACKAGE"))
  expect_no_condition(abort_packages_not_installed("rlang"))
})
