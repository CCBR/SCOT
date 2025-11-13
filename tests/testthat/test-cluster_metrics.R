test_that("clustering metric is consistent for (BRCA)", {
  
  brca.data = selectData("wu_et_al_BRCA")
  
  # run clustering
  brca.clustered = seurat_clustering(so_in = brca.data,
                                     npcs_in = 30)
  
  clusters = c("SCT_snn_res.0.8")
  
  res_sil = cluster_metrics(so = brca.clustered,
                                   cluster_list = clusters,
                                   dims = 1:20,
                                   reduction = "pca",
                                   silhouette = TRUE)
  
  # check proper dimnames and columns
  expect_true(is.matrix(res_sil))
  expect_equal(rownames(res_sil), c("clust_good", "clust_perm"))
  expect_equal(colnames(res_sil), c("CalinskiHarabasz", "DaviesBouldin", "Silhouette"))
  
  # no NAs present and matrix is numeric
  expect_true(all(is.finite(res_sil)))
  expect_true(is.numeric(res_sil))
})