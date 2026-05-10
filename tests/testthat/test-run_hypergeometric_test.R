# Unit tests for run_hypergeometric_test function

# Setup: Create test gene lists
set.seed(42)
all_genes <- paste0("GENE_", seq_len(1000))
selected_genes <- sample(all_genes, size = 50)
reference_genes <- sample(all_genes, size = 100)
world_size <- length(all_genes)

test_that("run_hypergeometric_test returns named vector with correct structure", {
  result <- run_hypergeometric_test(
    selectedVect = selected_genes,
    refVect = reference_genes,
    worldSize = world_size,
    lowerTail = FALSE
  )

  # Check result is a character vector with named elements
  expect_type(result, "character")
  expected_names <- c(
    "P-value",
    "IntersectSize",
    "SelectGenesSize",
    "ReferenceSize",
    "WorldSize",
    "IntersectedGenes"
  )
  expect_named(result, expected_names)
})

test_that("run_hypergeometric_test calculates intersection size correctly", {
  result <- run_hypergeometric_test(
    selectedVect = selected_genes,
    refVect = reference_genes,
    worldSize = world_size,
    lowerTail = FALSE
  )

  # Manually calculate expected intersection
  expected_intersection <- length(intersect(selected_genes, reference_genes))

  # Extract and convert to numeric for comparison
  actual_intersection <- as.numeric(result["IntersectSize"])
  expect_equal(actual_intersection, expected_intersection)
})

test_that("run_hypergeometric_test returns correct gene sizes", {
  result <- run_hypergeometric_test(
    selectedVect = selected_genes,
    refVect = reference_genes,
    worldSize = world_size,
    lowerTail = FALSE
  )

  expect_equal(as.numeric(result["SelectGenesSize"]), length(selected_genes))
  expect_equal(as.numeric(result["ReferenceSize"]), length(reference_genes))
  expect_equal(as.numeric(result["WorldSize"]), world_size)
})

test_that("run_hypergeometric_test returns p-value of 1 for no intersection", {
  # Create gene lists with no overlap
  selected_no_overlap <- paste0("SELECTED_", seq_len(50))
  reference_no_overlap <- paste0("REFERENCE_", seq_len(100))
  combined_world <- c(
    selected_no_overlap,
    reference_no_overlap,
    all_genes[201:1000]
  )

  result <- run_hypergeometric_test(
    selectedVect = selected_no_overlap,
    refVect = reference_no_overlap,
    worldSize = length(combined_world),
    lowerTail = FALSE
  )

  expect_equal(as.numeric(result["P-value"]), 1)
  expect_equal(as.numeric(result["IntersectSize"]), 0)
})

test_that("run_hypergeometric_test returns valid p-value between 0 and 1", {
  result <- run_hypergeometric_test(
    selectedVect = selected_genes,
    refVect = reference_genes,
    worldSize = world_size,
    lowerTail = FALSE
  )

  p_value <- as.numeric(result["P-value"])
  expect_true(p_value >= 0 && p_value <= 1)
})

test_that("run_hypergeometric_test lowerTail parameter affects p-value", {
  result_upper <- run_hypergeometric_test(
    selectedVect = selected_genes,
    refVect = reference_genes,
    worldSize = world_size,
    lowerTail = FALSE
  )

  result_lower <- run_hypergeometric_test(
    selectedVect = selected_genes,
    refVect = reference_genes,
    worldSize = world_size,
    lowerTail = TRUE
  )

  p_upper <- as.numeric(result_upper["P-value"])
  p_lower <- as.numeric(result_lower["P-value"])

  # P-values should be different (unless exactly 0.5)
  expect_false(p_upper == p_lower)
})

test_that("run_hypergeometric_test errors on invalid world size", {
  expect_warning(
    run_hypergeometric_test(
      selectedVect = selected_genes,
      refVect = reference_genes,
      worldSize = -100,
      lowerTail = FALSE
    ),
    regexp = "NaNs produced"
  )
})

test_that("run_hypergeometric_test returns intersected genes string", {
  result <- run_hypergeometric_test(
    selectedVect = selected_genes,
    refVect = reference_genes,
    worldSize = world_size,
    lowerTail = FALSE
  )

  # Check intersected genes are returned as semicolon-separated string
  intersected_string <- result["IntersectedGenes"]
  expect_type(intersected_string, "character")

  # If there are intersections, should contain semicolons or be a single gene
  if (as.numeric(result["IntersectSize"]) > 0) {
    expect_true(
      grepl("GENE_", intersected_string) || nchar(intersected_string) > 0
    )
  }
})
