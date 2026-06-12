# Tests for split_featurePlot

brca_data <- selectData("wu_et_al_BRCA")
features <- rownames(brca_data@assays$SCT@scale.data)[1:2]

test_that("split_featurePlot returns arranged ggplots in list format", {
  res <- split_featurePlot(
    so = brca_data,
    features = features,
    split_ident = "orig.ident",
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
  expect_true(all(vapply(res, function(x) inherits(x, "ggplot"), logical(1))))
})

test_that("split_featurePlot computes ncol/nrow when NA", {
  res <- split_featurePlot(
    so = brca_data,
    features = features,
    split_ident = "orig.ident",
    ncol = NA,
    nrow = NA,
    min.cutoff = "q5",
    max.cutoff = "q95",
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )

  # for multiple genes, check that a list of ggplots are returned
  expect_true(is.list(res))
  expect_true(all(vapply(res, function(x) inherits(x, "ggplot"), logical(1))))
})

test_that("split_featurePlot runs when provided with only nrow", {
  so <- brca_data
  features <- rownames(so@assays$SCT@scale.data)[1]
  res <- split_featurePlot(
    so = so,
    features = features,
    split_ident = "orig.ident",
    nrow = 1,
    ncol = NA,
    min.cutoff = "q1",
    max.cutoff = "q99",
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )
  expect_true(is.list(res))
  expect_true(inherits(res[[features]], "ggplot"))
})

test_that("split_featurePlot errors on invalid inputs", {
  expect_error(split_featurePlot(
    so,
    features = rownames(so)[1],
    split_ident = "bad_col"
  ))

  # Missing feature
  expect_error(split_featurePlot(
    so,
    features = "MissingGene",
    split_ident = "orig.ident"
  ))
})

test_that("split_featurePlot returns NULL silently with return_list = FALSE", {
  so <- brca_data
  res <- split_featurePlot(
    so = so,
    features = features[1],
    split_ident = "orig.ident",
    plot_image = FALSE,
    return_list = FALSE,
    reduction = "umap"
  )
  expect_null(res)
})

test_that("split_featurePlot runs when provided with only ncol", {
  so <- brca_data
  single_feature <- rownames(so@assays$SCT@scale.data)[1]

  res <- split_featurePlot(
    so = so,
    features = single_feature,
    split_ident = "orig.ident",
    ncol = 2,
    nrow = NA,
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )

  expect_true(is.list(res))
  expect_true(inherits(res[[single_feature]], "ggplot"))
})

test_that("split_featurePlot works with scalar numeric cutoff values", {
  so <- brca_data
  single_feature <- rownames(so@assays$SCT@scale.data)[1]

  res <- split_featurePlot(
    so = so,
    features = single_feature,
    split_ident = "orig.ident",
    min.cutoff = 0,
    max.cutoff = 5,
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )

  expect_true(is.list(res))
  expect_true(inherits(res[[single_feature]], "ggplot"))
})

test_that("split_featurePlot runs with label = TRUE", {
  so <- brca_data
  single_feature <- rownames(so@assays$SCT@scale.data)[1]

  res <- split_featurePlot(
    so = so,
    features = single_feature,
    split_ident = "orig.ident",
    label = TRUE,
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )

  expect_true(is.list(res))
  expect_true(inherits(res[[single_feature]], "ggplot"))
})

test_that("split_featurePlot works with order = TRUE", {
  so <- brca_data
  single_feature <- rownames(so@assays$SCT@scale.data)[1]

  res <- split_featurePlot(
    so = so,
    features = single_feature,
    split_ident = "orig.ident",
    order = TRUE,
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )

  expect_true(is.list(res))
  expect_true(inherits(res[[single_feature]], "ggplot"))
})

test_that("split_featurePlot works with slot = 'data'", {
  so <- brca_data
  single_feature <- rownames(so@assays$SCT@data)[1]

  res <- split_featurePlot(
    so = so,
    features = single_feature,
    split_ident = "orig.ident",
    slot = "data",
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )

  expect_true(is.list(res))
  expect_true(inherits(res[[single_feature]], "ggplot"))
})

test_that("split_featurePlot works with single split group", {
  so <- brca_data
  single_feature <- rownames(so@assays$SCT@scale.data)[1]

  # Create metadata with single group value
  so$single_group <- "group_1"

  res <- split_featurePlot(
    so = so,
    features = single_feature,
    split_ident = "single_group",
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )

  expect_true(is.list(res))
  expect_true(inherits(res[[single_feature]], "ggplot"))
})


test_that("split_featurePlot works with numeric split_ident values", {
  so <- brca_data
  single_feature <- rownames(so@assays$SCT@scale.data)[1]

  # Create numeric metadata
  so$numeric_group <- sample(1:3, ncol(so), replace = TRUE)

  res <- split_featurePlot(
    so = so,
    features = single_feature,
    split_ident = "numeric_group",
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )

  expect_true(is.list(res))
  expect_true(inherits(res[[single_feature]], "ggplot"))
})

test_that("split_featurePlot handles empty feature list", {
  so <- brca_data

  res <- split_featurePlot(
    so = so,
    features = character(0),
    split_ident = "orig.ident",
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )

  expect_true(is.list(res))
  expect_equal(length(res), 0)
})

test_that("split_featurePlot with asymmetric cutoff pair", {
  so <- brca_data
  single_feature <- rownames(so@assays$SCT@scale.data)[1]

  res <- split_featurePlot(
    so = so,
    features = single_feature,
    split_ident = "orig.ident",
    min.cutoff = "q10",
    max.cutoff = 3.5,
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )

  expect_true(is.list(res))
  expect_true(inherits(res[[single_feature]], "ggplot"))
})
