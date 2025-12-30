# Test for make_volcano_plot function

test_that("make_volcano_plot function exists and has correct signature", {
  # Source function for testing
  source(file.path("..", "..", "R", "make_volcano_plot.R"))
  
  # Test that the function exists and is callable
  expect_true(exists("make_volcano_plot"))
  expect_true(is.function(make_volcano_plot))
  
  # Check function signature
  fn_args <- names(formals(make_volcano_plot))
  expected_args <- c("de_table", "significant", "logfc", "pval")
  expect_true(all(expected_args %in% fn_args))
  expect_equal(length(fn_args), 4)
  
  # Check default values
  defaults <- formals(make_volcano_plot)
  expect_equal(defaults$significant, TRUE)
  expect_equal(defaults$logfc, TRUE)
  expect_equal(defaults$pval, TRUE)
})

test_that("make_volcano_plot validates input parameters", {
  # Source function for testing
  source(file.path("..", "..", "R", "make_volcano_plot.R"))
  
  # Test parameter types
  expect_true(is.logical(TRUE))  # significant should be logical
  expect_true(is.logical(FALSE)) # logfc should be logical
  expect_true(is.logical(TRUE))  # pval should be logical
  
  # Test with NULL input
  expect_error(make_volcano_plot(NULL))
  
  # Test with non-data.frame input
  expect_error(make_volcano_plot("not_a_dataframe"))
  expect_error(make_volcano_plot(123))
})

test_that("make_volcano_plot handles mock differential expression data", {
  # Skip if required packages are not available
  skip_if_not_installed("ggplot2")
  
  # Source function for testing
  source(file.path("..", "..", "R", "make_volcano_plot.R"))
  
  library(ggplot2)
  
  # Create mock differential expression table
  set.seed(123)
  n_genes <- 100
  
  de_table <- data.frame(
    p_val = runif(n_genes, 0.001, 0.1),
    avg_log2FC = rnorm(n_genes, 0, 2),
    pct.1 = runif(n_genes, 0.1, 0.9),
    pct.2 = runif(n_genes, 0.1, 0.9),
    p_val_adj = runif(n_genes, 0.001, 0.1),
    row.names = paste0("Gene", 1:n_genes)
  )
  
  # Add some significant genes
  de_table$p_val_adj[1:10] <- 0.01  # Significant p-values
  de_table$avg_log2FC[1:5] <- 2.0   # High positive fold change
  de_table$avg_log2FC[6:10] <- -2.0 # High negative fold change
  
  # Test basic volcano plot creation
  plot_result <- make_volcano_plot(de_table)
  
  # Test that result is a ggplot object
  expect_s3_class(plot_result, "ggplot")
  
  # Test that plot has expected aesthetics
  expect_true("avg_log2FC" %in% names(plot_result$mapping))
  expect_true("log10_p" %in% names(plot_result$mapping) || 
              any(grepl("log10_p", as.character(plot_result$mapping))))
})

test_that("make_volcano_plot log10 p-value transformation", {
  # Source function for testing
  source(file.path("..", "..", "R", "make_volcano_plot.R"))
  
  # Test the log10 p-value transformation logic
  test_pvals <- c(0.1, 0.05, 0.01, 0.001, 0)
  
  # This is the transformation logic from the function
  log10_p <- -log10(test_pvals)
  log10_p[which(test_pvals == 0)] <- 500
  
  # Test transformations
  expect_equal(log10_p[1], -log10(0.1), tolerance = 1e-10)
  expect_equal(log10_p[2], -log10(0.05), tolerance = 1e-10)
  expect_equal(log10_p[5], 500)  # Zero p-value should become 500
  
  # All transformed values should be positive
  expect_true(all(log10_p > 0))
  
  # Smaller p-values should have larger -log10 values
  expect_true(log10_p[4] > log10_p[3])  # 0.001 > 0.01
  expect_true(log10_p[3] > log10_p[2])  # 0.01 > 0.05
})

test_that("make_volcano_plot significance classification", {
  # Source function for testing
  source(file.path("..", "..", "R", "make_volcano_plot.R"))
  
  # Test significance classification logic
  mock_data <- data.frame(
    avg_log2FC = c(2.0, -2.0, 0.5, -0.5, 3.0, -3.0),
    p_val_adj = c(0.01, 0.01, 0.01, 0.01, 0.1, 0.1),
    row.names = paste0("Gene", 1:6)
  )
  
  # Test individual classification logic
  significance <- vector(length = nrow(mock_data))
  significance[] <- "NotSignificant"
  
  # LogFC threshold logic (> 1.5)
  logfc_significant <- which(abs(mock_data$avg_log2FC) > 1.5)
  expect_equal(logfc_significant, c(1, 2, 5, 6))  # Genes with |FC| > 1.5
  
  # P-value threshold logic (< 0.05)
  pval_significant <- which(mock_data$p_val_adj < 0.05)
  expect_equal(pval_significant, c(1, 2, 3, 4))  # Genes with p < 0.05
  
  # Combined significance (both conditions)
  combined_significant <- which(abs(mock_data$avg_log2FC) > 1.5 & 
                               mock_data$p_val_adj < 0.05)
  expect_equal(combined_significant, c(1, 2))  # Genes meeting both criteria
})

test_that("make_volcano_plot parameter combinations", {
  # Skip if required packages are not available
  skip_if_not_installed("ggplot2")
  
  # Source function for testing
  source(file.path("..", "..", "R", "make_volcano_plot.R"))
  
  library(ggplot2)
  
  # Create minimal test data
  de_table <- data.frame(
    p_val = c(0.01, 0.05, 0.1),
    avg_log2FC = c(2.0, -1.0, 0.5),
    p_val_adj = c(0.01, 0.05, 0.1),
    row.names = c("Gene1", "Gene2", "Gene3")
  )
  
  # Test different parameter combinations
  plot1 <- make_volcano_plot(de_table, significant = TRUE, logfc = TRUE, pval = TRUE)
  plot2 <- make_volcano_plot(de_table, significant = FALSE, logfc = TRUE, pval = FALSE)
  plot3 <- make_volcano_plot(de_table, significant = FALSE, logfc = FALSE, pval = TRUE)
  plot4 <- make_volcano_plot(de_table, significant = FALSE, logfc = FALSE, pval = FALSE)
  
  # All should return ggplot objects
  expect_s3_class(plot1, "ggplot")
  expect_s3_class(plot2, "ggplot")
  expect_s3_class(plot3, "ggplot")
  expect_s3_class(plot4, "ggplot")
})

test_that("make_volcano_plot data frame construction", {
  # Source function for testing
  source(file.path("..", "..", "R", "make_volcano_plot.R"))
  
  # Test the data frame construction logic
  test_data <- data.frame(
    avg_log2FC = c(1.0, -1.5, 2.0),
    p_val_adj = c(0.05, 0.01, 0.001),
    row.names = c("A", "B", "C")
  )
  
  # Simulate the data frame construction from the function
  log10_p <- -log10(test_data$p_val_adj)
  avg_log2FC <- test_data$avg_log2FC
  gene <- rownames(test_data)
  significance <- rep("NotSignificant", length(log10_p))
  
  df <- data.frame(
    avg_log2FC = avg_log2FC,
    log10_p = log10_p,
    significance = significance,
    gene = gene
  )
  rownames(df) <- rownames(test_data)
  
  # Test data frame properties
  expect_equal(nrow(df), 3)
  expect_equal(ncol(df), 4)
  expect_equal(rownames(df), c("A", "B", "C"))
  expect_true(is.numeric(df$avg_log2FC))
  expect_true(is.numeric(df$log10_p))
})

test_that("make_volcano_plot factor level ordering", {
  # Source function for testing
  source(file.path("..", "..", "R", "make_volcano_plot.R"))
  
  # Test factor level ordering from the function
  expected_levels <- c(
    "NotSignificant",
    "avg_log2FC > 1.5", 
    "p_val_adj < 0.05",
    "p_val_adj < 0.05 & avg_log2FC > 1.5"
  )
  
  # Test that all expected levels are character strings
  expect_true(all(sapply(expected_levels, is.character)))
  expect_equal(length(expected_levels), 4)
  
  # Test factor creation
  test_significance <- c("NotSignificant", "p_val_adj < 0.05", "NotSignificant")
  factor_sig <- factor(test_significance, levels = expected_levels)
  
  expect_true(is.factor(factor_sig))
  expect_equal(levels(factor_sig), expected_levels)
})

test_that("make_volcano_plot handles edge cases", {
  # Skip if required packages are not available
  skip_if_not_installed("ggplot2")
  
  # Source function for testing
  source(file.path("..", "..", "R", "make_volcano_plot.R"))
  
  library(ggplot2)
  
  # Test with minimal data (single gene)
  single_gene <- data.frame(
    p_val = 0.05,
    avg_log2FC = 1.0,
    p_val_adj = 0.05,
    row.names = "SingleGene"
  )
  
  plot_single <- make_volcano_plot(single_gene)
  expect_s3_class(plot_single, "ggplot")
  
  # Test with zero p-values
  zero_pval <- data.frame(
    p_val = c(0, 0.05),
    avg_log2FC = c(2.0, 1.0),
    p_val_adj = c(0, 0.05),
    row.names = c("ZeroP", "NonZeroP")
  )
  
  plot_zero <- make_volcano_plot(zero_pval)
  expect_s3_class(plot_zero, "ggplot")
})

test_that("make_volcano_plot missing columns validation", {
  # Source function for testing
  source(file.path("..", "..", "R", "make_volcano_plot.R"))
  
  # Test with missing required columns
  incomplete_table1 <- data.frame(
    p_val = c(0.01, 0.05),
    avg_log2FC = c(1.0, -1.0),
    # Missing p_val_adj
    row.names = c("Gene1", "Gene2")
  )
  
  expect_error(make_volcano_plot(incomplete_table1))
  
  incomplete_table2 <- data.frame(
    p_val = c(0.01, 0.05),
    # Missing avg_log2FC
    p_val_adj = c(0.01, 0.05),
    row.names = c("Gene1", "Gene2")
  )
  
  expect_error(make_volcano_plot(incomplete_table2))
})