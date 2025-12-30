context("run_singleR_db (targeted)")

test_that("run_singleR_db has expected signature", {
  expect_true(exists("run_singleR_db"))
  fmls <- formals(run_singleR_db)
  expect_setequal(names(fmls), c("so_in", "species"))
})

test_that("unknown species leaves object unchanged (meta.data)", {
  skip_if_not_installed("ontoProc")
  skip_if_not_installed("Seurat")

  suppressPackageStartupMessages(library(Seurat))

  set.seed(1)
  mat <- matrix(rpois(50, lambda = 3), nrow = 10, ncol = 5,
                dimnames = list(paste0("Gene", seq_len(10)), paste0("Cell", seq_len(5))))
  so <- CreateSeuratObject(counts = mat)
  md_before <- colnames(so@meta.data)

  so2 <- run_singleR_db(so, species = "unknown")
  expect_s4_class(so2, "Seurat")
  expect_identical(colnames(so2@meta.data), md_before)
})

test_that("human branch contains double-call bug for BlueprintEncodeData", {
  fn <- get("run_singleR_db")
  src <- paste(deparse(body(fn)), collapse = "\n")
  expect_match(src, "(celldex::BlueprintEncodeData())()", fixed = TRUE)
})
