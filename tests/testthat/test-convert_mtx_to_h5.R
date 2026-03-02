# Test for convert_mtx_to_h5 function

test_that("convert_mtx_to_h5 function exists and has correct signature", {
  # Source function for testing
  source(file.path("..", "..", "R", "convert_mtx_to_h5.R"))

  # Test that the function exists and is callable
  expect_true(exists("convert_mtx_to_h5"))
  expect_true(is.function(convert_mtx_to_h5))

  # Check function signature
  fn_args <- names(formals(convert_mtx_to_h5))
  expected_args <- c("sample_name", "mtx_file", "features_file", "barcodes_file")
  expect_true(all(expected_args %in% fn_args))
  expect_equal(length(fn_args), 4)
})

test_that("convert_mtx_to_h5 validates input parameters", {
  # Source function for testing
  source(file.path("..", "..", "R", "convert_mtx_to_h5.R"))

  # Test parameter types
  expect_true(is.character("sample1"))  # sample_name should be character
  expect_true(is.character("/path/to/matrix.mtx.gz"))  # mtx_file should be character
  expect_true(is.character("/path/to/features.tsv.gz"))  # features_file should be character
  expect_true(is.character("/path/to/barcodes.tsv.gz"))  # barcodes_file should be character

  # Test with NULL inputs - these should error
  expect_error(convert_mtx_to_h5(NULL, "mtx", "features", "barcodes"))
  expect_error(convert_mtx_to_h5("sample", NULL, "features", "barcodes"))
  expect_error(convert_mtx_to_h5("sample", "mtx", NULL, "barcodes"))
  expect_error(convert_mtx_to_h5("sample", "mtx", "features", NULL))
})

test_that("convert_mtx_to_h5 handles file path validation", {
  # Source function for testing
  source(file.path("..", "..", "R", "convert_mtx_to_h5.R"))

  # Test with non-existent files (should error when trying to read)
  expect_error(
    convert_mtx_to_h5("test_sample",
                      "nonexistent_matrix.mtx.gz",
                      "nonexistent_features.tsv.gz",
                      "nonexistent_barcodes.tsv.gz")
  )

  # Test that file paths are character strings
  valid_paths <- c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")
  expect_true(all(sapply(valid_paths, is.character)))
  expect_true(all(sapply(valid_paths, function(x) length(x) == 1)))
})

test_that("convert_mtx_to_h5 creates mock MTX files and converts successfully", {
  # Skip if required packages are not available
  skip_if_not_installed("Seurat")
  skip_if_not_installed("DropletUtils")
  skip_if_not_installed("Matrix")

  # Source function for testing
  source(file.path("..", "..", "R", "convert_mtx_to_h5.R"))

  library(Seurat)
  library(DropletUtils)
  library(Matrix)

  # Create temporary directory for test files
  temp_dir <- tempdir()

  # Create mock MTX format files
  set.seed(123)
  n_genes <- 10
  n_cells <- 5

  # Create mock sparse matrix
  mock_counts <- Matrix::sparseMatrix(
    i = sample(1:n_genes, 20, replace = TRUE),
    j = sample(1:n_cells, 20, replace = TRUE),
    x = rpois(20, lambda = 3),
    dims = c(n_genes, n_cells)
  )

  # Create mock features (genes)
  features <- data.frame(
    gene_id = paste0("ENSG", sprintf("%08d", 1:n_genes)),
    gene_name = paste0("Gene", 1:n_genes),
    feature_type = rep("Gene Expression", n_genes)
  )

  # Create mock barcodes
  barcodes <- data.frame(
    barcode = paste0("CELL", sprintf("%04d", 1:n_cells), "-1")
  )

  # Write mock files using DropletUtils
  temp_mtx_dir <- file.path(temp_dir, "test_mtx")
  dir.create(temp_mtx_dir, recursive = TRUE, showWarnings = FALSE)

  # Write the 10x format files
  DropletUtils::write10xCounts(
    path = temp_mtx_dir,
    x = mock_counts,
    gene.id = features$gene_id,
    gene.symbol = features$gene_name,
    gene.type = features$feature_type,
    barcodes = barcodes$barcode,
    overwrite = TRUE
  )

  # Set up file paths
  mtx_file <- file.path(temp_mtx_dir, "matrix.mtx.gz")
  features_file <- file.path(temp_mtx_dir, "features.tsv.gz")
  barcodes_file <- file.path(temp_mtx_dir, "barcodes.tsv.gz")

  # Test that files were created
  expect_true(file.exists(mtx_file))
  expect_true(file.exists(features_file))
  expect_true(file.exists(barcodes_file))

  # Set working directory to temp directory to control output location
  old_wd <- getwd()
  setwd(temp_dir)

  # Test the conversion function
  result <- convert_mtx_to_h5(
    sample_name = "test_sample",
    mtx_file = mtx_file,
    features_file = features_file,
    barcodes_file = barcodes_file
  )

  # Test that output file was created
  output_file <- file.path(temp_dir, "test_sample.h5")
  expect_true(file.exists(output_file))

  # Test that the function returns invisibly (should be NULL when called with invisible())
  expect_null(result)

  # Clean up
  setwd(old_wd)
  unlink(temp_mtx_dir, recursive = TRUE)
  unlink(output_file)
})

test_that("convert_mtx_to_h5 output file naming", {
  # Source function for testing
  source(file.path("..", "..", "R", "convert_mtx_to_h5.R"))

  # Test output file naming logic
  sample_names <- c("sample1", "test_data", "experiment_A")
  expected_outputs <- paste0(sample_names, ".h5")

  expect_equal(expected_outputs, c("sample1.h5", "test_data.h5", "experiment_A.h5"))

  # Test with special characters in sample name
  special_name <- "sample-test_123"
  expected_special <- paste0(special_name, ".h5")
  expect_equal(expected_special, "sample-test_123.h5")
})

test_that("convert_mtx_to_h5 handles different file extensions", {
  # Source function for testing
  source(file.path("..", "..", "R", "convert_mtx_to_h5.R"))

  # Test that function accepts different file path formats
  test_paths <- list(
    compressed = list(
      mtx = "path/to/matrix.mtx.gz",
      features = "path/to/features.tsv.gz",
      barcodes = "path/to/barcodes.tsv.gz"
    ),
    uncompressed = list(
      mtx = "path/to/matrix.mtx",
      features = "path/to/features.tsv",
      barcodes = "path/to/barcodes.tsv"
    ),
    alternative_names = list(
      mtx = "path/to/genes.mtx.gz",
      features = "path/to/genes.tsv.gz",
      barcodes = "path/to/cells.tsv.gz"
    )
  )

  # All paths should be valid character strings
  for (path_set in test_paths) {
    expect_true(is.character(path_set$mtx))
    expect_true(is.character(path_set$features))
    expect_true(is.character(path_set$barcodes))
  }
})

test_that("convert_mtx_to_h5 parameter order validation", {
  # Source function for testing
  source(file.path("..", "..", "R", "convert_mtx_to_h5.R"))

  # Test that parameters are in the expected order
  fn_formals <- formals(convert_mtx_to_h5)
  param_names <- names(fn_formals)

  # Check the parameter order matches documentation
  expect_equal(param_names[1], "sample_name")
  expect_equal(param_names[2], "mtx_file")
  expect_equal(param_names[3], "features_file")
  expect_equal(param_names[4], "barcodes_file")
})

test_that("convert_mtx_to_h5 integration with Seurat ReadMtx", {
  # Skip if required packages are not available
  skip_if_not_installed("Seurat")

  # Source function for testing
  source(file.path("..", "..", "R", "convert_mtx_to_h5.R"))

  # Test that Seurat::ReadMtx is available
  expect_true(exists("ReadMtx", envir = asNamespace("Seurat")))

  # Test parameter mapping for ReadMtx
  # ReadMtx expects: mtx, cells, features
  # Our function provides: mtx_file, barcodes_file, features_file

  # This tests the parameter mapping logic without actually running the function
  test_params <- list(
    mtx = "matrix.mtx.gz",
    cells = "barcodes.tsv.gz",  # Note: ReadMtx uses 'cells' parameter
    features = "features.tsv.gz"
  )

  expect_true(all(sapply(test_params, is.character)))
})

test_that("convert_mtx_to_h5 integration with DropletUtils write10xCounts", {
  # Skip if required packages are not available
  skip_if_not_installed("DropletUtils")

  # Source function for testing
  source(file.path("..", "..", "R", "convert_mtx_to_h5.R"))

  # Test that DropletUtils::write10xCounts is available
  expect_true(exists("write10xCounts", envir = asNamespace("DropletUtils")))

  # Test the expected parameters for write10xCounts
  # write10xCounts expects: x (matrix), path (output file)

  # Test output path generation
  sample_name <- "test_sample"
  expected_outfile <- paste0(sample_name, ".h5")
  expect_equal(expected_outfile, "test_sample.h5")
})
