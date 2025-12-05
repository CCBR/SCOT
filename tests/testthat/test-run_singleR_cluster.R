test_that("object_size", {
  pbmc <- suppressMessages(import_pbmc())
  set.seed(42)
  expect_equal(ncol(pbmc), 2638)
})

test_that("cluster_count", {
  pbmc <- suppressMessages(import_pbmc())
  set.seed(42)
  expect_equal(length(unique(pbmc$seurat_clusters)), 12)
})

test_that("cluster_cell_annot_type_count", {
  pbmc <- suppressMessages(import_pbmc())
  ref <- fetch_celldex_ref("hpca")
  Seurat::Idents(pbmc) <- "seurat_clusters"
  set.seed(42)
  singleR_cluster <- run_singleR_cluster(
    so_in = pbmc,
    ref_file = ref,
    label = "label.main"
  )
  expect_equal(length(unique(singleR_cluster)), 5)
  expect_equal(length(which(singleR_cluster == "B_cell")), 345, tolerance = 5)
  expect_equal(length(which(singleR_cluster == "Monocyte")), 672, tolerance = 5)
  expect_equal(length(which(singleR_cluster == "NK_cell")), 157, tolerance = 5)
  expect_equal(length(which(singleR_cluster == "T_cells")), 1452, tolerance = 5)
})
