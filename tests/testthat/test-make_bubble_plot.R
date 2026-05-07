brca_data <- load_fixture_data("wu_et_al_BRCA")
brca_genes <- head(rownames(brca_data), 5)
so_bubble <- brca_data

# set up custom ident for bubble plot testing
so_bubble$unit_test_ident <- rep(c("A", "B"), length.out = ncol(so_bubble))

bubble_plot <- make_bubble_plot(
  so = so_bubble,
  features = brca_genes,
  ident = "unit_test_ident"
)

test_that("bubble plot is successful for (BRCA)", {
  expect_true(inherits(bubble_plot, "ggplot"))
  expect_true(is.data.frame(bubble_plot$data))
  expect_true(all(
    c("Gene", "Group", "AvgExp", "PctExp") %in% colnames(bubble_plot$data)
  ))
})

test_that("bubble plot includes all feature and identity combinations", {
  bubble_plot <- make_bubble_plot(
    so = so_bubble,
    features = brca_genes,
    ident = "unit_test_ident"
  )

  n_groups <- length(unique(as.character(so_bubble$unit_test_ident)))
  expect_equal(nrow(bubble_plot$data), length(brca_genes) * n_groups)
  expect_true(all(
    bubble_plot$data$PctExp >= 0 & bubble_plot$data$PctExp <= 100
  ))
})

test_that("bubble plot supports scaled average expression", {
  bubble_plot_scaled <- make_bubble_plot(
    so = so_bubble,
    features = brca_genes,
    ident = "unit_test_ident",
    scale = TRUE
  )

  expect_true(inherits(bubble_plot_scaled, "ggplot"))
  expect_true(any(is.finite(bubble_plot_scaled$data$AvgExp)))
})

test_that("bubble plot respects custom ident columns", {
  so_custom <- brca_data
  so_custom$unit_test_ident <- rep(c("A", "B"), length.out = ncol(so_custom))

  bubble_plot_custom <- make_bubble_plot(
    so = so_custom,
    features = brca_genes,
    ident = "unit_test_ident"
  )

  expect_setequal(
    unique(as.character(bubble_plot_custom$data$Group)),
    c("A", "B")
  )
  expect_equal(nrow(bubble_plot_custom$data), length(brca_genes) * 2)
})

test_that("bubble plot errors on invalid ident", {
  expect_error(
    make_bubble_plot(
      so = brca_data,
      features = brca_genes,
      ident = "not_a_real_ident"
    )
  )
})

test_that("bubble plot errors when requested features are absent", {
  expect_error(
    make_bubble_plot(
      so = so_bubble,
      features = c("NOT_A_REAL_GENE_1", "NOT_A_REAL_GENE_2"),
      ident = "unit_test_ident"
    )
  )
})
