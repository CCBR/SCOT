brca_data <- load_fixture_data("wu_et_al_BRCA")

test_that("filter_doublets returns a Seurat object", {
  filtered_result <- filter_doublets(
    brca_data,
    doublet_finder_method = "scDblFinder"
  )
  expect_s4_class(filtered_result, "Seurat")
})

test_that("unknown method errors clearly", {
  expect_error(
    filter_doublets(brca_data, doublet_finder_method = "bogus"),
    regexp = "No valid method selected"
  )
})

test_that("scDblFinder branch annotates and subsets singlets (if available)", {
  set.seed(7)
  test_counts <- matrix(
    rpois(400, lambda = 3),
    nrow = 20,
    ncol = 20,
    dimnames = list(paste0("Gene", 1:20), paste0("Cell", 1:20))
  )
  test_seurat <- Seurat::CreateSeuratObject(counts = test_counts)

  filtered_result <- suppressWarnings(filter_doublets(
    test_seurat,
    doublet_finder_method = "scDblFinder"
  ))
  expect_s4_class(filtered_result, "Seurat")
  expect_true("scDblFinder_label" %in% colnames(filtered_result@meta.data))
  expect_lte(ncol(filtered_result), ncol(test_seurat))
})
