# Tests for split_featurePlot

brca_data <- load_fixture_data("wu_et_al_BRCA")
top_var_features <- c("S100P", "CCL3", "KRT17", "KRT19", "FCN3")

test_that("split_featurePlot returns arranged ggplots in list format", {
  res <- split_featurePlot(
    so = brca_data,
    features = top_var_features,
    split_ident = "Phase",
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
  expect_equal(sort(names(res)), sort(top_var_features))
  expect_true(all(vapply(res, function(x) inherits(x, "ggplot"), logical(1))))

  # # Check that plots have data and layers
  # for (feature_plot in res) {
  #   expect_true(nrow(feature_plot$data) > 0)
  #   expect_true(length(feature_plot$layers) > 0)
  # }
})

test_that("split_featurePlot computes ncol/nrow when NA", {
  res <- split_featurePlot(
    so = brca_data,
    features = top_var_features,
    split_ident = "Phase",
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

  # Check that the plots have faceted structure (split by split_ident)
  for (feature_plot in res) {
    expect_true(!is.null(feature_plot$facet))
  }
})

test_that("split_featurePlot runs when provided with only nrow", {
  so <- brca_data
  features <- top_var_features[1]
  res <- split_featurePlot(
    so = so,
    features = features,
    split_ident = "Phase",
    nrow = 1,
    ncol = NA,
    min.cutoff = "q1",
    max.cutoff = "q99",
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )
  expect_true(is.list(res))

  # Access the plot directly by feature name
  plot <- res[[features]]
  expect_true(inherits(plot, "ggplot"))

  # Check plot structure - don't access $data directly
  expect_true(length(plot$layers) >= 1)
  expect_true(!is.null(plot$facet))
  expect_true(!is.null(plot$mapping))
})

test_that("split_featurePlot errors on invalid inputs", {
  # Invalid split_ident column
  expect_error(
    split_featurePlot(
      brca_data,
      features = "S100P",
      split_ident = "bad_col"
    ),
    "not found|bad_col"
  )

  # Missing feature
  expect_error(
    split_featurePlot(
      brca_data,
      features = "MissingGene",
      split_ident = "Phase"
    ),
    "MissingGene|not found|feature"
  )
})

test_that("split_featurePlot returns NULL silently with return_list = FALSE", {
  so <- brca_data
  res <- split_featurePlot(
    so = so,
    features = top_var_features[1],
    split_ident = "Phase",
    plot_image = FALSE,
    return_list = FALSE,
    reduction = "umap"
  )
  expect_null(res)
})

test_that("split_featurePlot runs when provided with only ncol", {
  so <- brca_data
  single_feature <- top_var_features[1]

  res <- split_featurePlot(
    so = so,
    features = single_feature,
    split_ident = "Phase",
    ncol = 2,
    nrow = NA,
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )

  expect_true(is.list(res))
  expect_true(inherits(res[[single_feature]], "ggplot"))

  # Check plot structure - don't access $data directly
  plot <- res[[single_feature]]
  expect_true(length(plot$layers) >= 1)
  expect_true(!is.null(plot$facet))
  expect_true(!is.null(plot$mapping))
})

test_that("split_featurePlot works with scalar numeric cutoff values", {
  so <- brca_data
  single_feature <- top_var_features[1]

  res <- split_featurePlot(
    so = so,
    features = single_feature,
    split_ident = "Phase",
    min.cutoff = 0,
    max.cutoff = 5,
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )

  expect_true(is.list(res))
  expect_true(inherits(res[[single_feature]], "ggplot"))

  # Verify plot structure without accessing $data directly
  plot <- res[[single_feature]]
  expect_true(length(plot$layers) >= 1)
  expect_true(!is.null(plot$facet))
  expect_true(!is.null(plot$mapping))
})

test_that("split_featurePlot runs with label = TRUE", {
  so <- brca_data
  single_feature <- top_var_features[1]

  res <- split_featurePlot(
    so = so,
    features = single_feature,
    split_ident = "Phase",
    label = TRUE,
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )

  expect_true(is.list(res))
  expect_true(inherits(res[[single_feature]], "ggplot"))

  # Verify plot structure without accessing $data directly
  plot <- res[[single_feature]]
  expect_true(length(plot$layers) >= 1)
  expect_true(!is.null(plot$facet))
  expect_true(!is.null(plot$mapping))
})

test_that("split_featurePlot works with order = TRUE", {
  so <- brca_data
  single_feature <- top_var_features[1]

  res <- split_featurePlot(
    so = so,
    features = single_feature,
    split_ident = "Phase",
    order = TRUE,
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )

  expect_true(is.list(res))
  expect_true(inherits(res[[single_feature]], "ggplot"))

  # Check plot structure
  plot <- res[[single_feature]]
  expect_true(length(plot$layers) >= 1)
  expect_true(!is.null(plot$facet))
  expect_true(!is.null(plot$mapping))
})

test_that("split_featurePlot works with slot = 'data'", {
  so <- brca_data
  single_feature <- top_var_features[1]

  res <- split_featurePlot(
    so = so,
    features = single_feature,
    split_ident = "Phase",
    slot = "data",
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )

  expect_true(is.list(res))
  expect_true(inherits(res[[single_feature]], "ggplot"))

  # Check plot has expected structure
  plot <- res[[single_feature]]
  expect_true(length(plot$layers) >= 1)
  expect_true(!is.null(plot$facet))
  expect_true(!is.null(plot$mapping))
})

test_that("split_featurePlot works with single split group", {
  so <- brca_data
  single_feature <- top_var_features[1]

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
  single_feature <- top_var_features[1]

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
    split_ident = "Phase",
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )

  expect_true(is.list(res))
  expect_equal(length(res), 0)
})

test_that("split_featurePlot with asymmetric cutoff pair", {
  so <- brca_data
  single_feature <- top_var_features[1]

  res <- split_featurePlot(
    so = so,
    features = single_feature,
    split_ident = "Phase",
    min.cutoff = "q10",
    max.cutoff = 3.5,
    plot_image = FALSE,
    return_list = TRUE,
    reduction = "umap"
  )

  expect_true(is.list(res))
  expect_true(inherits(res[[single_feature]], "ggplot"))

  # Verify plot structure with mixed cutoff types
  plot <- res[[single_feature]]
  expect_true(length(plot$layers) >= 1)
  expect_true(!is.null(plot$facet))
  expect_true(!is.null(plot$mapping))
})
