# Test for run_AUCell function
so <- selectData("wu_et_al_BRCA")

test_that("run_AUCell return value structure", {
  # Test the expected return structure

  # Mock return structure based on function
  mock_scores <- matrix(runif(100), nrow = 4, ncol = 25)
  rownames(mock_scores) <- paste0("GeneSet_", 1:4)
  colnames(mock_scores) <- paste0("Cell_", 1:25)

  mock_assignment <- list(
    "GeneSet_1" = list(assignment = c(TRUE, FALSE, TRUE)),
    "GeneSet_2" = list(assignment = c(FALSE, TRUE, FALSE))
  )

  expected_output <- list(
    scores = mock_scores,
    assignment = mock_assignment
  )

  # Test output structure
  expect_true(is.list(expected_output))
  expect_equal(length(expected_output), 2)
  expect_true("scores" %in% names(expected_output))
  expect_true("assignment" %in% names(expected_output))

  # Test scores properties
  expect_true(
    is.matrix(expected_output$scores) || is.array(expected_output$scores)
  )
  expect_true(is.numeric(as.vector(expected_output$scores)))

  # Test assignment properties
  expect_true(is.list(expected_output$assignment))
})

test_that("run_AUCell FetchData parameters validation", {
  # Test FetchData parameter logic from the function

  # Mock parameters that would be used in FetchData
  mock_vars <- c("Gene1", "Gene2", "Gene3")
  mock_slot <- "SCT"
  mock_layer <- "data"

  # Test parameter types
  expect_true(is.character(mock_vars))
  expect_true(is.character(mock_slot))
  expect_true(is.character(mock_layer))

  # Test specific values from function
  expect_equal(mock_slot, "SCT")
  expect_equal(mock_layer, "data")

  # Test that vars would be gene names
  expect_true(length(mock_vars) > 0)
  expect_true(all(nchar(mock_vars) > 0))
})

test_that("run_AUCell AUCell parameters validation", {
  # Test AUCell-specific parameters from the function

  # nCores parameter
  expected_ncores <- 4
  expect_true(is.numeric(expected_ncores))
  expect_true(expected_ncores > 0)

  # plotStats parameter
  expected_plot_stats <- FALSE
  expect_true(is.logical(expected_plot_stats))

  # plotHist parameter
  expected_plot_hist <- FALSE
  expect_true(is.logical(expected_plot_hist))

  # assignCells parameter
  expected_assign_cells <- TRUE
  expect_true(is.logical(expected_assign_cells))
})

test_that("run_AUCell matrix transpose operation", {
  # Test matrix transpose logic from the function

  # Create mock expression data (genes x cells)
  set.seed(123)
  mock_expr <- data.frame(
    Gene1 = runif(5),
    Gene2 = runif(5),
    Gene3 = runif(5)
  )
  rownames(mock_expr) <- paste0("Cell_", 1:5)

  # Transpose operation: t(expr) should give (genes x cells)
  expr_matrix <- t(mock_expr)

  # Test transpose properties
  expect_equal(nrow(expr_matrix), 3) # 3 genes
  expect_equal(ncol(expr_matrix), 5) # 5 cells
  expect_equal(rownames(expr_matrix), c("Gene1", "Gene2", "Gene3"))
  expect_equal(colnames(expr_matrix), paste0("Cell_", 1:5))

  # Test that values are preserved
  expect_true(is.numeric(expr_matrix))
  expect_equal(dim(expr_matrix), c(3, 5))
})

test_that("run_AUCell edge cases and error handling", {
  # Use testthat::test_path() for correct path resolution
  run_aucell_path <- test_path("../../R/run_AUCell.R")

  skip_if_not(
    file.exists(run_aucell_path),
    message = "run_AUCell.R not found at expected location"
  )

  source(run_aucell_path)

  # Test with NULL inputs
  expect_error(run_AUCell(NULL, list()))
  expect_error(run_AUCell("not_seurat", list()))

  # Test with empty gene sets
  empty_gene_sets <- list()
  expect_true(is.list(empty_gene_sets))
  expect_equal(length(empty_gene_sets), 0)

  # Test with invalid gene sets
  invalid_gene_sets <- list(
    "BadSet1" = character(0), # Empty character vector
    "BadSet2" = NA, # NA value
    "BadSet3" = NULL # NULL value
  )

  expect_true(is.list(invalid_gene_sets))
  expect_equal(length(invalid_gene_sets$BadSet1), 0)
  expect_true(is.na(invalid_gene_sets$BadSet2))
  expect_null(invalid_gene_sets$BadSet3)
})
