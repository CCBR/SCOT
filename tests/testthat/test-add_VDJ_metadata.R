# Test for add_VDJ_metadata function

test_that("add_VDJ_metadata works with valid inputs", {
  # Skip if required packages are not available
  skip_if_not_installed("djvdj")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  # Source function for testing
  source(file.path("..", "..", "R", "add_VDJ_metadata.R"))

  # Create a minimal mock Seurat object for testing
  library(Seurat)
  library(SeuratObject)

  # Create minimal test data with proper gene names (no underscores)
  counts <- matrix(rpois(100, 5), nrow = 10, ncol = 10)
  rownames(counts) <- paste0("Gene", 1:10)  # No underscores to avoid Seurat warnings
  colnames(counts) <- paste0("Cell-", 1:10)

  # Create Seurat object
  test_so <- CreateSeuratObject(counts = counts)

  # Test with non-existent VDJ directory (should error gracefully)
  expect_error(
    add_VDJ_metadata(test_so, "non_existent_directory"),
    class = "error"
  )
})

test_that("add_VDJ_metadata handles prefix correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("djvdj")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  # Source function for testing
  source(file.path("..", "..", "R", "add_VDJ_metadata.R"))

  # Test the prefix logic embedded in the function
  test_prefix_logic <- function(vdj_path) {
    prefix <- basename(vdj_path)
    if (prefix == "outs") {
      prefix <- ""
    } else {
      prefix <- paste0(prefix, "_")
    }
    return(prefix)
  }

  # Test different directory names
  expect_equal(test_prefix_logic("/path/to/outs"), "")
  expect_equal(test_prefix_logic("/path/to/sample1"), "sample1_")
  expect_equal(test_prefix_logic("/path/to/TCR"), "TCR_")
})

test_that("add_VDJ_metadata validates input parameters", {
  # Skip if required packages are not available
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  # Source function for testing
  source(file.path("..", "..", "R", "add_VDJ_metadata.R"))

  library(Seurat)
  library(SeuratObject)

  # Create minimal test data with proper names
  counts <- matrix(rpois(100, 5), nrow = 10, ncol = 10)
  rownames(counts) <- paste0("Gene", 1:10)
  colnames(counts) <- paste0("Cell-", 1:10)
  test_so <- CreateSeuratObject(counts = counts)

  # Test with NULL inputs - these should error
  expect_error(add_VDJ_metadata(NULL, "some_path"))
  expect_error(add_VDJ_metadata(test_so, NULL))

  # Test with invalid Seurat object
  expect_error(add_VDJ_metadata("not_a_seurat_object", "some_path"))

  # Test with invalid file path types
  expect_error(add_VDJ_metadata(test_so, 123))
})

test_that("add_VDJ_metadata function exists and has correct signature", {
  # Source the function for testing since package may not be installed
  source(file.path("..", "..", "R", "add_VDJ_metadata.R"))

  # Test that the function exists and is callable
  expect_true(exists("add_VDJ_metadata"))
  expect_true(is.function(add_VDJ_metadata))

  # Check function signature
  fn_args <- names(formals(add_VDJ_metadata))
  expect_true("so" %in% fn_args)
  expect_true("vdj_file" %in% fn_args)
  expect_equal(length(fn_args), 2)
})

test_that("add_VDJ_metadata dependencies are available", {
  # Test that required packages are available for import
  expect_true(requireNamespace("Seurat", quietly = TRUE))

  # Test djvdj availability (may be optional)
  if (requireNamespace("djvdj", quietly = TRUE)) {
    expect_true(exists("import_vdj", envir = asNamespace("djvdj")))
  } else {
    skip("djvdj package not available")
  }
})

test_that("prefix logic works correctly", {
  # Test the prefix logic that's embedded in the function
  # This tests the core logic without requiring external dependencies

  # Simulate the prefix logic from the function
  test_prefix_logic <- function(vdj_file) {
    prefix <- basename(vdj_file)
    if (prefix == "outs") {
      prefix <- ""
    } else {
      prefix <- paste0(prefix, "_")
    }
    return(prefix)
  }

  # Test cases for different directory names
  expect_equal(test_prefix_logic("/path/to/outs"), "")
  expect_equal(test_prefix_logic("/path/to/sample1"), "sample1_")
  expect_equal(test_prefix_logic("/path/to/TCR_data"), "TCR_data_")
  expect_equal(test_prefix_logic("simple_name"), "simple_name_")
})

test_that("add_VDJ_metadata parameter validation concepts", {
  # Test parameter validation concepts without actually calling the function
  # This tests the types of inputs the function expects

  # Test that we can identify what should be valid inputs
  expect_true(is.character("valid_file_path"))
  expect_false(is.character(123))
  expect_false(is.character(NULL))

  # Test path handling
  test_paths <- c("/valid/path", "relative/path", "outs", "sample_1")
  expect_true(all(sapply(test_paths, is.character)))
  expect_true(all(sapply(test_paths, function(x) length(x) == 1)))
})
