# move inputs to top of file to avoid redundant code in test_that blocks
brca_data <- load_fixture_data("wu_et_al_BRCA")
input_cells <- colnames(brca_data)
use_npcs <- 30

test_that("seurat_clustering produces expected output structure (BRCA)", {
  set.seed(42)
  brca_clustered = seurat_clustering(so_in = brca_data, npcs_in = use_npcs)

  expect_s4_class(brca_clustered, "Seurat")

  # check that pca, umap, and cluster information is present
  reduction_names <- names(brca_clustered@reductions)
  expect_true("pca" %in% tolower(reduction_names))
  expect_true(ncol(Seurat::Embeddings(brca_clustered[["pca"]])) %in% c(1:50))

  expect_true("umap" %in% tolower(reduction_names))
  expect_equal(ncol(Seurat::Embeddings(brca_clustered[["umap"]])), 3)

  expect_true("seurat_clusters" %in% colnames(brca_clustered@meta.data))
  expect_true(is.factor(SeuratObject::Idents(brca_clustered)))
  expect_equal(
    length(SeuratObject::Idents(brca_clustered)),
    ncol(brca_clustered)
  )
  expect_identical(colnames(brca_clustered), input_cells)

  # Nearest-neighbor graph exists
  expect_true(length(brca_clustered@graphs) >= 1)

  set.seed(1)
  low_res <- seurat_clustering(
    brca_clustered,
    npcs_in = 10,
    resolution = 0.2,
    algorithm = 1
  )
  set.seed(2)
  high_res <- seurat_clustering(
    brca_clustered,
    npcs_in = 10,
    resolution = 1.2,
    algorithm = 1
  )

  # check number of high resolution clusters >= number of high resolution clusters
  ncl_low <- nlevels(SeuratObject::Idents(low_res))
  ncl_high <- nlevels(SeuratObject::Idents(high_res))

  expect_true(ncl_high >= ncl_low)
})

test_that("invalid inputs error clearly", {
  expect_error(
    seurat_clustering(brca_data, npcs_in = 0),
    regexp = "npcs_in"
  )
  expect_error(
    seurat_clustering(brca_data, npcs_in = 2.5),
    regexp = "npcs_in"
  )
  expect_error(
    seurat_clustering(brca_data, npcs_in = 5, resolution = 0),
    regexp = "resolution"
  )
  expect_error(
    seurat_clustering(brca_data, npcs_in = 5, algorithm = 4),
    regexp = "algorithm"
  )
})
