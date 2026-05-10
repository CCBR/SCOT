ifnb <- readRDS(testthat::test_path("fixtures/ifnb_sub.rds"))

preprocessed_result <- preprocess_sample(
  so_in = ifnb,
  species = "hg38",
  npcs_in = 10,
  num_ctrl = 100,
  replace = TRUE
)

test_that("preprocess_sample returns a Seurat object with human species", {
  expect_s4_class(preprocessed_result, "Seurat")
  expect_true(ncol(preprocessed_result) > 0)
  expect_true(nrow(preprocessed_result) > 0)
})

test_that("preprocess_sample adds cell cycle scoring metadata", {
  # Check that cell cycle scoring columns are added
  metadata_cols <- colnames(preprocessed_result@meta.data)
  expect_true("S.Score" %in% metadata_cols)
  expect_true("G2M.Score" %in% metadata_cols)
  expect_true("Phase" %in% metadata_cols)
})

test_that("preprocess_sample adds dimensionality reductions", {
  # Check that PCA and UMAP reductions are added
  reductions <- names(preprocessed_result@reductions)
  expect_true("pca" %in% reductions)
  expect_true("umap" %in% reductions)
})

test_that("preprocess_sample adds clustering results", {
  # Check that clustering columns are added
  metadata_cols <- colnames(preprocessed_result@meta.data)
  cluster_cols <- metadata_cols[grepl("res\\.", metadata_cols)]
  expect_true(length(cluster_cols) > 0)
})

test_that("preprocess_sample errors on invalid species", {
  expect_error(
    preprocess_sample(
      so_in = ifnb,
      species = "invalid_sp",
      npcs_in = 10
    ),
    regexp = "species|hg19|hg38|mm10"
  )
})

test_that("preprocess_sample errors on non-Seurat object", {
  expect_error(
    preprocess_sample(
      so_in = "not_a_seurat",
      species = "hg38",
      npcs_in = 10
    ),
    regexp = "Seurat|object|class"
  )
})

test_that("preprocess_sample errors on NULL object", {
  expect_error(
    preprocess_sample(
      so_in = NULL,
      species = "hg38",
      npcs_in = 10
    ),
    regexp = "NULL|Seurat|object"
  )
})
