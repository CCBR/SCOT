brca.data <- selectData("wu_et_al_BRCA")
brca.clustered <- seurat_clustering(so_in = brca.data, npcs_in = 30)

test_that("clustering metric is consistent for (BRCA)", {
  clusters <- c("SCT_snn_res.0.2", "SCT_snn_res.0.8")

  res_sil <- cluster_metrics(
    so = brca.clustered,
    cluster_list = clusters,
    dims = 1:20,
    reduction = "pca",
    silhouette = TRUE
  )

  expect_true(is.matrix(res_sil))
  expect_equal(rownames(res_sil), clusters)
  expect_equal(
    colnames(res_sil),
    c("CalinskiHarabasz", "DaviesBouldin", "Silhouette")
  )

  # explicit silhouette inclusion check
  expect_true("Silhouette" %in% colnames(res_sil))

  expect_true(all(is.finite(res_sil)))
  expect_true(is.numeric(res_sil))
})

test_that("cluster_metrics returns 2 columns when silhouette is FALSE", {
  res <- cluster_metrics(
    so = brca.clustered,
    cluster_list = c("SCT_snn_res.0.2", "SCT_snn_res.0.8"),
    dims = 1:20,
    reduction = "pca",
    silhouette = FALSE
  )

  expect_true(is.matrix(res))
  expect_equal(colnames(res), c("CalinskiHarabasz", "DaviesBouldin"))
  expect_equal(ncol(res), 2)

  # explicit silhouette omission check
  expect_false("Silhouette" %in% colnames(res))
})

test_that("cluster label coercion works for numeric, character, and factor metadata columns", {
  so_types <- brca.clustered
  base_label <- "SCT_snn_res.0.2"
  base_clusters <- as.numeric(unlist(so_types[[base_label]]))

  so_types[["cluster_numeric_test"]] <- base_clusters
  so_types[["cluster_character_test"]] <- as.character(base_clusters)
  so_types[["cluster_factor_test"]] <- factor(base_clusters)

  cluster_labels <- c(
    "cluster_numeric_test",
    "cluster_character_test",
    "cluster_factor_test"
  )

  res_types <- cluster_metrics(
    so = so_types,
    cluster_list = cluster_labels,
    dims = 1:20,
    reduction = "pca",
    silhouette = FALSE
  )

  expect_true(is.matrix(res_types))
  expect_equal(rownames(res_types), cluster_labels)
  expect_equal(colnames(res_types), c("CalinskiHarabasz", "DaviesBouldin"))
  expect_true(all(is.finite(res_types)))

  # all encodings represent the same partition, so scores should match
  expect_equal(
    res_types["cluster_numeric_test", ],
    res_types["cluster_character_test", ],
    tolerance = 1e-8
  )
  expect_equal(
    res_types["cluster_numeric_test", ],
    res_types["cluster_factor_test", ],
    tolerance = 1e-8
  )
})

test_that("cluster_metrics errors on missing cluster column", {
  expect_error(
    cluster_metrics(
      so = brca.clustered,
      cluster_list = c("SCT_snn_res.0.2", "not_a_real_cluster"),
      dims = 1:20,
      reduction = "pca",
      silhouette = TRUE
    )
  )
})

test_that("cluster_metrics errors on invalid reduction", {
  expect_error(
    cluster_metrics(
      so = brca.clustered,
      cluster_list = c("SCT_snn_res.0.2"),
      dims = 1:20,
      reduction = "not_a_reduction",
      silhouette = TRUE
    )
  )
})

test_that("cluster_metrics errors when dims exceed available PCs", {
  expect_error(
    cluster_metrics(
      so = brca.clustered,
      cluster_list = c("SCT_snn_res.0.2"),
      dims = 1:200,
      reduction = "pca",
      silhouette = TRUE
    )
  )
})
