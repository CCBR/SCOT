context("filter_doublets")

test_that("filter_doublets has expected signature and default", {
  expect_true(exists("filter_doublets"))
  fmls <- formals(filter_doublets)
  expect_setequal(names(fmls), c("so_in", "doublet_finder_method"))
  expect_identical(fmls$doublet_finder_method, "DoubletFinder")
})

test_that("unknown method errors clearly", {
  skip_if_not_installed("Seurat")
  suppressPackageStartupMessages(library(Seurat))

  counts <- matrix(rpois(60, lambda = 2), nrow = 12, ncol = 5,
                   dimnames = list(paste0("Gene", 1:12), paste0("Cell", 1:5)))
  so <- CreateSeuratObject(counts = counts)
  expect_error(filter_doublets(so, doublet_finder_method = "bogus"), 
               regexp = "No valid method selected")
})

test_that("scDblFinder branch annotates and subsets singlets (if available)", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("scDblFinder")

  suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
  })

  set.seed(7)
  counts <- matrix(rpois(400, lambda = 3), nrow = 20, ncol = 20,
                   dimnames = list(paste0("Gene", 1:20), paste0("Cell", 1:20)))
  so <- CreateSeuratObject(counts = counts)

  res <- suppressWarnings(filter_doublets(so, doublet_finder_method = "scDblFinder"))
  expect_s4_class(res, "Seurat")
  expect_true("scDblFinder_label" %in% colnames(res@meta.data))
  expect_lte(ncol(res), ncol(so))
})

test_that("source contains DoubletFinder calls (static check)", {
  fn <- get("filter_doublets")
  src <- paste(deparse(body(fn)), collapse = "\n")
  expect_match(src, "DoubletFinder::doubletFinder(", fixed = TRUE)
  expect_match(src, "DoubletFinder::modelHomotypic(", fixed = TRUE)
  expect_match(src, "doubletFinder_label", fixed = TRUE)
})
