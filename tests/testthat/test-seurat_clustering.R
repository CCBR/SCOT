test_that("seurat_clustering has expected signature and defaults", {
  expect_true(exists("seurat_clustering"))
  fmls <- formals(seurat_clustering)
  expect_setequal(names(fmls), c("so_in", "npcs_in", "resolution", "algorithm"))
  expect_identical(fmls$resolution, 0.8)
  expect_identical(fmls$algorithm, 3)
})

test_that("clustering metric is consistent for (BRCA)", {
  brca.data = load_fixture_data("wu_et_al_BRCA")
  input_cells <- colnames(brca.data)

  set.seed(42)
  brca.clustered = seurat_clustering(so_in = brca.data, npcs_in = 30)

  expect_s4_class(brca.clustered, "Seurat")

  # check that pca, umap, and cluster information is present
  reduction_names <- names(brca.clustered@reductions)
  expect_true("pca" %in% tolower(reduction_names))
  expect_true(ncol(Seurat::Embeddings(brca.clustered[["pca"]])) %in% c(1:50))

  expect_true("umap" %in% tolower(reduction_names))
  expect_equal(ncol(Seurat::Embeddings(brca.clustered[["umap"]])), 3)

  expect_true("seurat_clusters" %in% colnames(brca.clustered@meta.data))
  expect_true(is.factor(SeuratObject::Idents(brca.clustered)))
  expect_equal(
    length(SeuratObject::Idents(brca.clustered)),
    ncol(brca.clustered)
  )
  expect_identical(colnames(brca.clustered), input_cells)

  # Nearest-neighbor graph exists
  expect_true(length(brca.clustered@graphs) >= 1)

  set.seed(1)
  low_res <- seurat_clustering(
    brca.clustered,
    npcs_in = 10,
    resolution = 0.2,
    algorithm = 1
  )
  set.seed(2)
  high_res <- seurat_clustering(
    brca.clustered,
    npcs_in = 10,
    resolution = 1.2,
    algorithm = 1
  )

  # check number of high resolution clusters >= number of high resolution clusters
  ncl_low <- nlevels(SeuratObject::Idents(low_res))
  ncl_high <- nlevels(SeuratObject::Idents(high_res))

  expect_true(ncl_high >= ncl_low)
})

test_that("resolution and algorithm arguments are forwarded to FindClusters", {
  fn <- get("seurat_clustering")
  src <- paste(deparse(body(fn)), collapse = "\n")
  expect_match(src, "resolution = resolution", fixed = TRUE)
  expect_match(src, "algorithm = algorithm", fixed = TRUE)
})

test_that("invalid inputs error clearly", {
  brca.data <- load_fixture_data("wu_et_al_BRCA")

  expect_error(
    seurat_clustering(brca.data, npcs_in = 0),
    regexp = "npcs_in"
  )
  expect_error(
    seurat_clustering(brca.data, npcs_in = 2.5),
    regexp = "npcs_in"
  )
  expect_error(
    seurat_clustering(brca.data, npcs_in = 5, resolution = 0),
    regexp = "resolution"
  )
  expect_error(
    seurat_clustering(brca.data, npcs_in = 5, algorithm = 4),
    regexp = "algorithm"
  )
})
