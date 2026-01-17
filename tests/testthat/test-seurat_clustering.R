test_that("clustering metric is consistent for (BRCA)", {

  brca.data = selectData("wu_et_al_BRCA")

  set.seed(42)
  brca.clustered = seurat_clustering(so_in = brca.data,
                                   npcs_in = 30)

  expect_s4_class(brca.clustered, "Seurat")

  # check that pca, umap, and cluster information is present
  expect_true("pca" %in% names(Reductions(brca.clustered)))
  expect_true(ncol(Embeddings(brca.clustered[["pca"]])) %in% c(1:50))

  expect_true("umap" %in% names(Reductions(brca.clustered)))
  expect_equal(ncol(Embeddings(brca.clustered[["umap"]])), 3)

  expect_true("seurat_clusters" %in% colnames(brca.clustered@meta.data))
  expect_true(is.factor(Idents(brca.clustered)))
  expect_equal(length(Idents(brca.clustered)), ncol(brca.clustered))

  # Nearest-neighbor graph exists
  expect_true(length(Graphs(brca.clustered)) >= 1)

  set.seed(1)
  low_res  <- seurat_clustering(brca.clustered, npcs_in = 10, resolution = 0.2, algorithm = 1)
  set.seed(2)
  high_res <- seurat_clustering(brca.clustered, npcs_in = 10, resolution = 1.2, algorithm = 1)

  # check number of high resolution clusters >= number of high resolution clusters
  ncl_low  <- nlevels(Idents(low_res))
  ncl_high <- nlevels(Idents(high_res))

  expect_true(ncl_high >= ncl_low)
})
