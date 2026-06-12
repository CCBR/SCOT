# Test for run_AUCell function
brca_data <- load_fixture_data("wu_et_al_BRCA")
available_genes <- rownames(brca_data)

test_gene_sets <- list(
  "test_set_1" = available_genes[1:3],
  "test_set_2" = available_genes[4:6]
)

aucell_result <- run_AUCell(
  so = brca_data,
  gene_sets = test_gene_sets
)

test_that("run_AUCell returns list with scores and assignment", {
  expect_type(aucell_result, "list")
  expect_true("scores" %in% names(aucell_result))
  expect_true("assignment" %in% names(aucell_result))
})

test_that("run_AUCell scores matrix has correct dimensions and values", {
  scores_matrix <- aucell_result$scores

  expect_true(is.matrix(scores_matrix))
  expect_equal(nrow(scores_matrix), length(test_gene_sets))
  expect_equal(ncol(scores_matrix), ncol(brca_data))
  expect_equal(rownames(scores_matrix), names(test_gene_sets))
  expect_true(is.numeric(scores_matrix))
  expect_true(all(is.finite(scores_matrix)))
  expect_true(all(scores_matrix >= 0 & scores_matrix <= 1))
})

test_that("run_AUCell assignment is list with expected structure", {
  assignment_list <- aucell_result$assignment

  expect_type(assignment_list, "list")
  expect_true("test_set_1" %in% names(assignment_list))
  expect_true("assignment" %in% names(assignment_list$test_set_1))
})

test_that("run_AUCell errors on invalid Seurat object", {
  test_gene_sets <- list("test_set" = available_genes[1:2])

  expect_error(
    run_AUCell(so = NULL, gene_sets = test_gene_sets),
    regexp = "Seurat|invalid|object"
  )
})

test_that("run_AUCell errors on empty gene sets", {
  expect_error(
    run_AUCell(so = brca_data, gene_sets = list()),
    regexp = "gene|set|empty"
  )
})
