# Test for make_volcano_plot function

test_that("make_volcano_plot handles mock differential expression data", {
  # Create mock differential expression table
  set.seed(123)
  n_genes <- 100

  de_table <- data.frame(
    p_val = runif(n_genes, 0.001, 0.1),
    avg_log_2_fc = rnorm(n_genes, 0, 2),
    pct.1 = runif(n_genes, 0.1, 0.9),
    pct.2 = runif(n_genes, 0.1, 0.9),
    p_val_adj = runif(n_genes, 0.001, 0.1),
    row.names = paste0("Gene", 1:n_genes)
  )

  # Add some significant genes
  de_table$p_val_adj[1:10] <- 0.01 # Significant p-values
  de_table$avg_log_2_fc[1:5] <- 2.0 # High positive fold change
  de_table$avg_log_2_fc[6:10] <- -2.0 # High negative fold change

  # Test basic volcano plot creation
  plot_result <- make_volcano_plot(de_table)

  # Test that result is a ggplot object
  expect_s3_class(plot_result, "ggplot")

  # Test that plot has the expected aesthetics (x, y, color)
  expect_true("x" %in% names(plot_result$mapping))
  expect_true("y" %in% names(plot_result$mapping))
  expect_true(
    "colour" %in%
      names(plot_result$mapping) ||
      "color" %in% names(plot_result$mapping)
  )

  # Test plot data shape and structure
  expect_true(is.data.frame(plot_result$data))
  expect_equal(nrow(plot_result$data), n_genes)
  expect_true(all(
    c("avg_log_2_fc", "log10_p", "significance", "gene") %in%
      colnames(plot_result$data)
  ))

  # Test that plot has layers (geom_point at minimum)
  expect_true(length(plot_result$layers) >= 1)
  expect_equal(class(plot_result$layers[[1]]$geom)[1], "GeomPoint")
})


test_that("make_volcano_plot parameter combinations", {
  # Create minimal test data
  de_table <- data.frame(
    p_val = c(0.01, 0.05, 0.1),
    avg_log_2_fc = c(2.0, -1.0, 0.5),
    p_val_adj = c(0.01, 0.05, 0.1),
    row.names = c("Gene1", "Gene2", "Gene3")
  )

  # Test different parameter combinations
  plot1 <- make_volcano_plot(
    de_table,
    significant = TRUE,
    logfc = TRUE,
    pval = TRUE
  )
  plot2 <- make_volcano_plot(
    de_table,
    significant = FALSE,
    logfc = TRUE,
    pval = FALSE
  )
  plot3 <- make_volcano_plot(
    de_table,
    significant = FALSE,
    logfc = FALSE,
    pval = TRUE
  )
  plot4 <- make_volcano_plot(
    de_table,
    significant = FALSE,
    logfc = FALSE,
    pval = FALSE
  )

  # All should return ggplot objects with proper structure
  for (plt in list(plot1, plot2, plot3, plot4)) {
    expect_s3_class(plt, "ggplot")
    expect_true(is.data.frame(plt$data))
    expect_equal(nrow(plt$data), 3)
    expect_true(all(
      c("avg_log_2_fc", "log10_p", "significance", "gene") %in%
        colnames(plt$data)
    ))
    expect_true(length(plt$layers) >= 1)
    expect_equal(class(plt$layers[[1]]$geom)[1], "GeomPoint")
  }
})


test_that("make_volcano_plot handles edge cases", {
  # Test with minimal data (single gene)
  single_gene <- data.frame(
    p_val = 0.05,
    avg_log_2_fc = 1.0,
    p_val_adj = 0.05,
    row.names = "SingleGene"
  )

  plot_single <- make_volcano_plot(single_gene)
  expect_s3_class(plot_single, "ggplot")
  expect_equal(nrow(plot_single$data), 1)
  expect_true(all(
    c("avg_log_2_fc", "log10_p", "significance", "gene") %in%
      colnames(plot_single$data)
  ))
  expect_true(length(plot_single$layers) >= 1)

  # Test with zero p-values
  zero_p_val <- data.frame(
    p_val = c(0, 0.05),
    avg_log_2_fc = c(2.0, 1.0),
    p_val_adj = c(0, 0.05),
    row.names = c("ZeroP", "NonZeroP")
  )

  plot_zero <- make_volcano_plot(zero_p_val)
  expect_s3_class(plot_zero, "ggplot")
  expect_equal(nrow(plot_zero$data), 2)
  expect_true(all(
    c("avg_log_2_fc", "log10_p", "significance", "gene") %in%
      colnames(plot_zero$data)
  ))
  expect_true(length(plot_zero$layers) >= 1)
})

test_that("make_volcano_plot missing columns validation", {
  # Test with missing required columns
  incomplete_table_1 <- data.frame(
    p_val = c(0.01, 0.05),
    avg_log_2_fc = c(1.0, -1.0),
    # Missing p_val_adj
    row.names = c("Gene1", "Gene2")
  )

  expect_error(
    make_volcano_plot(incomplete_table_1),
    "non-numeric|argument to mathematical function"
  )

  incomplete_table_2 <- data.frame(
    p_val = c(0.01, 0.05),
    # Missing avg_log_2_fc
    p_val_adj = c(0.01, 0.05),
    row.names = c("Gene1", "Gene2")
  )

  expect_error(
    make_volcano_plot(incomplete_table_2),
    "non-numeric|argument to mathematical function"
  )
})
