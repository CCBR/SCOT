test_that("convert_human_genes", {
  input_human <- c("TP53", "GAPDH", "ACTB")
  expected_mouse <- c("Trp53", "Gapdh", "Actb")

  expect_equal(convert_human_gene_list(input_human), expected_mouse)
})

test_that("human_unique_genes", {
  input_human <- c("TP53", "SRGAP2C", "HYDIN2", "ARGAP11B")
  expect_equal(convert_human_gene_list(input_human), "Trp53")
  #Error in mouse-side retrieval if ONLY human-unique genes
  expect_error(convert_human_gene_list(input_human[-1]))
})

test_that("no_valid_genes", {
  input_non_existent <- c("Aaaa", "abababa", "XKCD")
  expect_error(convert_human_gene_list(input_non_existent))
})
