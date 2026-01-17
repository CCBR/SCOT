# Tests for split_featurePlot

context("split_featurePlot")

test_that("split_featurePlot exists and has correct signature", {
  source(file.path("..", "..", "R", "split_featurePlot.R"))
  expect_true(exists("split_featurePlot"))
  expect_true(is.function(split_featurePlot))

  args <- names(formals(split_featurePlot))
  expected <- c(
    "so", "features", "split_ident", "label", "ncol", "nrow",
    "min.cutoff", "max.cutoff", "plot_image", "return_list",
    "slot", "order", "reduction"
  )
  expect_equal(args, expected)

  defaults <- formals(split_featurePlot)
  expect_equal(defaults$label, FALSE)
  expect_true(is.na(defaults$ncol))
  expect_true(is.na(defaults$nrow))
  expect_true(is.na(defaults$min.cutoff))
  expect_true(is.na(defaults$max.cutoff))
  expect_equal(defaults$plot_image, FALSE)
  expect_equal(defaults$return_list, FALSE)
  expect_equal(defaults$slot, "scale.data")
  expect_equal(defaults$order, FALSE)
  expect_null(defaults$reduction)
})

test_that("split_featurePlot dependencies available", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpubr")
  skip_if_not_installed("Matrix")
})

make_mock_seurat <- function(n_genes = 50, n_cells = 40) {
  library(Seurat)
  library(Matrix)
  set.seed(123)
  counts <- Matrix::sparseMatrix(
    i = sample(seq_len(n_genes), n_genes * 4, replace = TRUE),
    j = sample(seq_len(n_cells), n_cells * 4, replace = TRUE),
    x = rpois(n_genes * 4, 2),
    dims = c(n_genes, n_cells)
  )
  rownames(counts) <- paste0("Gene", seq_len(n_genes))
  colnames(counts) <- paste0("Cell", seq_len(n_cells))
  so <- Seurat::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
  so$group <- rep(c("A", "B"), length.out = n_cells)
  Seurat::DefaultAssay(so) <- "RNA"
  so <- Seurat::NormalizeData(so)
  so <- Seurat::ScaleData(so)
  so <- Seurat::FindVariableFeatures(so)
  so <- Seurat::RunPCA(so, npcs = 10)
  so <- Seurat::RunUMAP(so, dims = 1:10)
  so
}

test_that("split_featurePlot returns arranged plots with return_list", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpubr")
  skip_if_not_installed("Matrix")

  source(file.path("..", "..", "R", "split_featurePlot.R"))

  so <- make_mock_seurat(40, 30)
  features <- rownames(so)[1:2]
  res <- split_featurePlot(
    so = so,
    features = features,
    split_ident = "group",
    label = TRUE,
    ncol = NA,
    nrow = NA,
    min.cutoff = "q10",
    max.cutoff = "q90",
    plot_image = FALSE,
    return_list = TRUE,
    slot = "scale.data",
    order = FALSE,
    reduction = "umap"
  )
  expect_true(is.list(res))
  expect_equal(sort(names(res)), sort(features))
  # ggpubr::annotate_figure returns an object of class "ggpubr"
  expect_true(all(vapply(res, function(x) inherits(x, "ggpubr"), logical(1))))
})

test_that("split_featurePlot computes ncol/nrow when NA", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpubr")
  skip_if_not_installed("Matrix")

  source(file.path("..", "..", "R", "split_featurePlot.R"))
  so <- make_mock_seurat(30, 10)
  features <- rownames(so)[1]
  res <- split_featurePlot(
    so = so,
    features = features,
    split_ident = "group",
    ncol = NA,
    nrow = NA,
    min.cutoff = "q5",
    max.cutoff = "q95",
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )
  expect_true(is.list(res))
  expect_true(inherits(res[[features]], "ggpubr"))
})

test_that("split_featurePlot respects provided nrow only", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpubr")
  skip_if_not_installed("Matrix")

  source(file.path("..", "..", "R", "split_featurePlot.R"))
  so <- make_mock_seurat(20, 8)
  features <- rownames(so)[1]
  res <- split_featurePlot(
    so = so,
    features = features,
    split_ident = "group",
    nrow = 1,
    ncol = NA,
    min.cutoff = "q1",
    max.cutoff = "q99",
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )
  expect_true(is.list(res))
  expect_true(inherits(res[[features]], "ggpubr"))
})

test_that("split_featurePlot errors on invalid inputs", {
  source(file.path("..", "..", "R", "split_featurePlot.R"))
  expect_error(split_featurePlot(NULL, features = "Gene1", split_ident = "group"))

  # Invalid split_ident
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  so <- make_mock_seurat(20, 8)
  expect_error(split_featurePlot(so, features = rownames(so)[1], split_ident = "bad_col"))

  # Missing feature
  expect_error(split_featurePlot(so, features = "MissingGene", split_ident = "group"))
})
