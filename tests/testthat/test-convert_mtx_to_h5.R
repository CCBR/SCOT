# Unit tests for convert_mtx_to_h5 function

# Load fixture files from fixtures directory (relative to test file location)
mtx_file <- file.path("fixtures", "test_matrix.mtx.gz")
barcodes_file <- file.path("fixtures", "test_barcodes.tsv.gz")
features_file <- file.path("fixtures", "test_features.tsv.gz")

# Create temporary directory for test outputs
temp_fixture_dir <- tempdir()

test_that("convert_mtx_to_h5 creates h5 output file", {
  output_sample_name <- tempfile(tmpdir = temp_fixture_dir, fileext = "")
  output_h5_file <- paste0(output_sample_name, ".h5")

  # Run conversion
  result <- convert_mtx_to_h5(
    sample_name = output_sample_name,
    mtx_file = mtx_file,
    features_file = features_file,
    barcodes_file = barcodes_file
  )

  # Check that h5 file was created
  expect_true(file.exists(output_h5_file))

  # Cleanup
  if (file.exists(output_h5_file)) {
    file.remove(output_h5_file)
  }
})

test_that("convert_mtx_to_h5 errors on invalid mtx file path", {
  expect_error(
    convert_mtx_to_h5(
      sample_name = "test_sample",
      mtx_file = "/nonexistent/path/matrix.mtx.gz",
      features_file = features_file,
      barcodes_file = barcodes_file
    ),
    regexp = "cannot find||expression"
  )
})

test_that("convert_mtx_to_h5 errors on invalid features file path", {
  expect_error(
    convert_mtx_to_h5(
      sample_name = "test_sample",
      mtx_file = mtx_file,
      features_file = "/nonexistent/path/features.tsv.gz",
      barcodes_file = barcodes_file
    ),
    regexp = "cannot find|feature list"
  )
})

test_that("convert_mtx_to_h5 errors on invalid barcodes file path", {
  expect_error(
    convert_mtx_to_h5(
      sample_name = "test_sample",
      mtx_file = mtx_file,
      features_file = features_file,
      barcodes_file = "/nonexistent/path/barcodes.tsv.gz"
    ),
    regexp = "cannot find|barcode list"
  )
})

test_that("convert_mtx_to_h5 generates h5 with correct naming convention", {
  sample_name <- "test_sample_output"
  output_sample_name <- file.path(temp_fixture_dir, sample_name)
  expected_h5_file <- paste0(output_sample_name, ".h5")

  convert_mtx_to_h5(
    sample_name = output_sample_name,
    mtx_file = mtx_file,
    features_file = features_file,
    barcodes_file = barcodes_file
  )

  # Verify h5 file has correct name format
  expect_true(grepl("\\.h5$", expected_h5_file))
  expect_true(file.exists(expected_h5_file))

  # Cleanup
  if (file.exists(expected_h5_file)) {
    file.remove(expected_h5_file)
  }
})
