pkg_test <- function(pkgs) {
  return(abort_packages_not_installed(pkgs))
}
test_that("multiplication works", {
  expect_error(pkg_test("NOT_A_PACKAGE"))
  expect_no_condition(pkg_test("rlang"))
})
