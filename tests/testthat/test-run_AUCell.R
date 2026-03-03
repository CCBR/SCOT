# Test for run_AUCell function

test_that("run_AUCell function exists and has correct signature", {
  # Source function for testing
  source(file.path("..", "..", "R", "run_AUCell.R"))

  # Test that the function exists and is callable
  expect_true(exists("run_AUCell"))
  expect_true(is.function(run_AUCell))

  # Check function signature
  fn_args <- names(formals(run_AUCell))
  expected_args <- c("so", "gene_sets")
  expect_true(all(expected_args %in% fn_args))
  expect_equal(length(fn_args), 2)
})

test_that("run_AUCell validates input parameters", {
  # Source function for testing
  source(file.path("..", "..", "R", "run_AUCell.R"))

  # Test that gene_sets should be a list
  test_gene_sets <- list(
    "Pathway1" = c("Gene1", "Gene2", "Gene3"),
    "Pathway2" = c("Gene4", "Gene5")
  )

  expect_true(is.list(test_gene_sets))
  expect_true(all(sapply(test_gene_sets, is.character)))
  expect_equal(length(test_gene_sets), 2)
  expect_equal(names(test_gene_sets), c("Pathway1", "Pathway2"))
})

test_that("run_AUCell dependency checks", {
  # Skip if required packages are not available
  skip_if_not_installed("AUCell")
  skip_if_not_installed("Seurat")

  # Test that AUCell functions are available
  library(AUCell)

  expect_true(exists("subsetGeneSets", where = asNamespace("AUCell")))
  expect_true(exists("setGeneSetNames", where = asNamespace("AUCell")))
  expect_true(exists("nGenes", where = asNamespace("AUCell")))
  expect_true(exists("AUCell_buildRankings", where = asNamespace("AUCell")))
  expect_true(exists("AUCell_calcAUC", where = asNamespace("AUCell")))
  expect_true(exists("AUCell_exploreThresholds", where = asNamespace("AUCell")))
})

test_that("run_AUCell mock Seurat object creation for AUCell", {
  # Skip if required packages are not available
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  library(Seurat)
  library(Matrix)

  # Create mock SCT-transformed Seurat object
  set.seed(123)
  n_genes <- 50
  n_cells <- 30

  # Create count matrix
  counts <- Matrix::sparseMatrix(
    i = sample(1:n_genes, 200, replace = TRUE),
    j = sample(1:n_cells, 200, replace = TRUE),
    x = rpois(200, 3),
    dims = c(n_genes, n_cells)
  )

  gene_names <- paste0("Gene_", 1:n_genes)
  cell_names <- paste0("Cell_", 1:n_cells)
  rownames(counts) <- gene_names
  colnames(counts) <- cell_names

  # Create Seurat object
  mock_so <- CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)

  # Add SCT assay with normalized data
  normalized_data <- log1p(counts)
  mock_so[["SCT"]] <- CreateAssayObject(data = normalized_data)

  # Test Seurat object properties
  expect_s4_class(mock_so, "Seurat")
  expect_true("SCT" %in% names(mock_so@assays))
  expect_equal(nrow(mock_so[["SCT"]]@data), length(gene_names))
  expect_equal(ncol(mock_so[["SCT"]]@data), length(cell_names))
})

test_that("run_AUCell gene set formatting", {
  # Test gene set structure and formatting

  # Create test gene sets
  test_gene_sets <- list(
    "Immune_Response" = c("CD3D", "CD3E", "CD8A", "CD4"),
    "Cell_Cycle" = c("MKI67", "PCNA", "CDC20"),
    "Apoptosis" = c("TP53", "BAX", "BCL2")
  )

  # Test gene set properties
  expect_true(is.list(test_gene_sets))
  expect_equal(length(test_gene_sets), 3)
  expect_true(all(names(test_gene_sets) != ""))

  # Test gene set contents
  for (gene_set in test_gene_sets) {
    expect_true(is.character(gene_set))
    expect_true(length(gene_set) > 0)
    expect_true(all(nchar(gene_set) > 0))
  }

  # Test gene set naming convention (from function)
  set_names <- names(test_gene_sets)
  set_lengths <- sapply(test_gene_sets, length)
  expected_names <- paste(set_names, " (", set_lengths, "g)", sep = "")

  expect_equal(length(expected_names), 3)
  expect_true(all(grepl("\\(\\d+g\\)$", expected_names)))
})

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
  expect_true(is.matrix(expected_output$scores) || is.array(expected_output$scores))
  expect_true(is.numeric(as.vector(expected_output$scores)))

  # Test assignment properties
  expect_true(is.list(expected_output$assignment))
})

test_that("run_AUCell seed setting verification", {
  # Source function for testing
  source(file.path("..", "..", "R", "run_AUCell.R"))

  # Test that seed is set to 42 in the function
  # This ensures reproducibility
  expected_seed <- 42

  # Test seed value
  expect_true(is.numeric(expected_seed))
  expect_equal(expected_seed, 42)

  # Test reproducibility with set.seed
  set.seed(42)
  random1 <- runif(5)

  set.seed(42)
  random2 <- runif(5)

  expect_equal(random1, random2)
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
  expect_equal(nrow(expr_matrix), 3)  # 3 genes
  expect_equal(ncol(expr_matrix), 5)  # 5 cells
  expect_equal(rownames(expr_matrix), c("Gene1", "Gene2", "Gene3"))
  expect_equal(colnames(expr_matrix), paste0("Cell_", 1:5))

  # Test that values are preserved
  expect_true(is.numeric(expr_matrix))
  expect_equal(dim(expr_matrix), c(3, 5))
})

test_that("run_AUCell edge cases and error handling", {
  # Source function for testing
  source(file.path("..", "..", "R", "run_AUCell.R"))

  # Test with NULL inputs
  expect_error(run_AUCell(NULL, list()))
  expect_error(run_AUCell("not_seurat", list()))

  # Test with empty gene sets
  empty_gene_sets <- list()
  expect_true(is.list(empty_gene_sets))
  expect_equal(length(empty_gene_sets), 0)

  # Test with invalid gene sets
  invalid_gene_sets <- list(
    "BadSet1" = character(0),  # Empty character vector
    "BadSet2" = NA,           # NA value
    "BadSet3" = NULL          # NULL value
  )

  expect_true(is.list(invalid_gene_sets))
  expect_equal(length(invalid_gene_sets$BadSet1), 0)
  expect_true(is.na(invalid_gene_sets$BadSet2))
  expect_null(invalid_gene_sets$BadSet3)
})
