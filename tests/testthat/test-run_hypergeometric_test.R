context("run_hypergeometric_test")

# Ensure function is available when running file directly
if (!exists("run_hypergeometric_test")) {
  suppressWarnings(
    source("/Users/bianjh/Documents/Projects/Development/SCOT/R/run_hypergeometric_test.R")
  )
}

test_that("run_hypergeometric_test has expected signature", {
  expect_true(exists("run_hypergeometric_test"))
  fmls <- formals(run_hypergeometric_test)
  expect_setequal(names(fmls), c("selectedVect", "refVect", "worldSize", "lowerTail"))
  expect_identical(fmls$lowerTail, FALSE)
})

test_that("returns named character vector with correct fields", {
  sel <- c("A", "B", "C", "D")
  ref <- c("B", "D", "X")
  world <- 100

  out <- run_hypergeometric_test(sel, ref, world)
  expect_type(out, "character")
  expect_setequal(names(out), c(
    "P-value", "IntersectSize", "SelectGenesSize",
    "ReferenceSize", "WorldSize", "IntersectedGenes"
  ))
  expect_identical(out[["IntersectSize"]], as.character(2L))
  expect_identical(out[["SelectGenesSize"]], as.character(length(sel)))
  expect_identical(out[["ReferenceSize"]], as.character(length(ref)))
  expect_identical(out[["WorldSize"]], as.character(world))
  expect_identical(out[["IntersectedGenes"]], "B;D")
})

test_that("p-value matches phyper computation (upper tail)", {
  sel <- c("a", "b", "c", "d", "e")
  ref <- c("b", "d", "f", "g")
  world <- 200
  inter <- length(intersect(sel, ref))
  pv_expected <- stats::phyper(q = inter, m = length(ref), n = world - length(ref),
                               k = length(sel), lower.tail = FALSE)

  out <- run_hypergeometric_test(sel, ref, world)
  expect_equal(as.numeric(out[["P-value"]]), pv_expected)
})

test_that("p-value equals 1 when intersection size is 0", {
  sel <- c("a", "b")
  ref <- c("x", "y")
  world <- 50
  out <- run_hypergeometric_test(sel, ref, world)
  expect_identical(out[["P-value"]], as.character(1))
  expect_identical(out[["IntersectedGenes"]], "")
})

test_that("lowerTail=TRUE flips tail behavior", {
  sel <- c("a", "b", "c", "d")
  ref <- c("b", "z")
  world <- 120
  inter <- length(intersect(sel, ref))
  pv_upper <- stats::phyper(q = inter, m = length(ref), n = world - length(ref),
                            k = length(sel), lower.tail = FALSE)
  pv_lower <- stats::phyper(q = inter, m = length(ref), n = world - length(ref),
                            k = length(sel), lower.tail = TRUE)

  out_upper <- run_hypergeometric_test(sel, ref, world, lowerTail = FALSE)
  out_lower <- run_hypergeometric_test(sel, ref, world, lowerTail = TRUE)
  expect_equal(as.numeric(out_upper[["P-value"]]), pv_upper)
  expect_equal(as.numeric(out_lower[["P-value"]]), pv_lower)
})
