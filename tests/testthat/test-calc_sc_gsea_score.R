# Test for calc_sc_gsea_score function

test_that("calc_sc_gsea_score computes scores correctly", {
  # Create a mock differential expression table
  deTable <- data.frame(
    p_val = c(0.001, 0.05, 0.1, 0.0001),
    avg_log2FC = c(2.5, -1.2, 0.8, -3.1),
    pct.1 = c(0.8, 0.3, 0.6, 0.9),
    pct.2 = c(0.2, 0.7, 0.4, 0.1),
    p_val_adj = c(0.01, 0.15, 0.25, 0.001),
    row.names = c("Gene1", "Gene2", "Gene3", "Gene4")
  )

  # Calculate scores
  scores <- calc_sc_gsea_score(deTable)

  # Test basic properties
  expect_type(scores, "double")
  expect_equal(length(scores), nrow(deTable))
  expect_equal(names(scores), rownames(deTable))

  # Test that scores are named correctly
  expect_true(all(names(scores) %in% rownames(deTable)))

  # Test formula components - manually calculate expected values for verification
  # For Gene1: sign(2.5) * -log10(0.001) * max(0.8, 0.2) + 2.5
  # = 1 * 3 * 0.8 + 2.5 = 4.9
  expected_gene1 <- 1 * (-log10(0.001)) * max(0.8, 0.2) + 2.5
  expect_equal(as.numeric(scores["Gene1"]), expected_gene1, tolerance = 1e-10)
})

test_that("calc_sc_gsea_score handles zero p-values correctly", {
  # Create table with zero p-values
  deTable <- data.frame(
    p_val = c(0, 0.05, 0),
    avg_log2FC = c(2.0, -1.5, -2.5),
    pct.1 = c(0.9, 0.4, 0.7),
    pct.2 = c(0.1, 0.6, 0.3),
    p_val_adj = c(0, 0.1, 0),
    row.names = c("ZeroGene1", "NonZeroGene", "ZeroGene2")
  )

  scores <- calc_sc_gsea_score(deTable)

  # For zero p-values, formula should be: sign(log2FC) * 500 * max(pct.1, pct.2) + log2FC
  # ZeroGene1: sign(2.0) * 500 * max(0.9, 0.1) + 2.0 = 1 * 500 * 0.9 + 2.0 = 452
  expected_zero1 <- 1 * 500 * max(0.9, 0.1) + 2.0
  expect_equal(
    as.numeric(scores["ZeroGene1"]),
    expected_zero1,
    tolerance = 1e-10
  )

  # ZeroGene2: sign(-2.5) * 500 * max(0.7, 0.3) + (-2.5) = -1 * 500 * 0.7 + (-2.5) = -352.5
  expected_zero2 <- -1 * 500 * max(0.7, 0.3) + (-2.5)
  expect_equal(
    as.numeric(scores["ZeroGene2"]),
    expected_zero2,
    tolerance = 1e-10
  )
})

test_that("calc_sc_gsea_score handles edge cases", {
  # Test with minimal data - single gene
  deTable_minimal <- data.frame(
    p_val = 0.5,
    avg_log2FC = 1.0,
    pct.1 = 0.5,
    pct.2 = 0.3,
    p_val_adj = 0.6,
    row.names = "SingleGene"
  )

  scores_minimal <- calc_sc_gsea_score(deTable_minimal)
  expect_length(scores_minimal, 1)
  expect_equal(names(scores_minimal), "SingleGene")

  # Test with zero fold change
  deTable_zero_fc <- data.frame(
    p_val = c(0.01, 0.05),
    avg_log2FC = c(0, 0),
    pct.1 = c(0.6, 0.4),
    pct.2 = c(0.4, 0.6),
    p_val_adj = c(0.05, 0.1),
    row.names = c("ZeroFC1", "ZeroFC2")
  )

  scores_zero_fc <- calc_sc_gsea_score(deTable_zero_fc)
  expect_length(scores_zero_fc, 2)
  # When log2FC is 0, sign should be 0, so the first part becomes 0, leaving only avg_log2FC
  expect_equal(as.numeric(scores_zero_fc["ZeroFC1"]), 0, tolerance = 1e-10)

  # Test with NA values - should result in NA scores
  deTable_na <- data.frame(
    p_val = c(0.01, NA),
    avg_log2FC = c(1.5, 2.0),
    pct.1 = c(0.7, 0.8),
    pct.2 = c(0.3, 0.2),
    p_val_adj = c(0.05, NA),
    row.names = c("Gene1", "NA_Gene2")
  )
  scores_na <- calc_sc_gsea_score(deTable_na)
  expect_true(!is.na(scores_na["Gene1"]))
  expect_true(is.na(scores_na["NA_Gene2"]))
})

test_that("calc_sc_gsea_score validates input structure", {
  # Test with missing required columns
  incomplete_table <- data.frame(
    p_val = c(0.01, 0.05),
    avg_log2FC = c(1.0, -0.5),
    # Missing pct.1 and pct.2
    row.names = c("Gene1", "Gene2")
  )

  expect_error(calc_sc_gsea_score(incomplete_table))

  # Test with NULL input
  expect_error(calc_sc_gsea_score(NULL))

  # Test with non-data.frame input
  expect_error(calc_sc_gsea_score("not_a_dataframe"))
})

test_that("calc_sc_gsea_score handles various p-value ranges", {
  # Test with very small and large p-values
  deTable <- data.frame(
    p_val = c(1e-10, 0.99, 1e-300, 0.5),
    avg_log2FC = c(3.0, -1.0, 2.5, -0.5),
    pct.1 = c(0.8, 0.2, 0.9, 0.4),
    pct.2 = c(0.1, 0.8, 0.1, 0.6),
    p_val_adj = c(1e-8, 0.99, 1e-250, 0.7),
    row.names = c("VerySmallP", "LargeP", "TinyP", "MediumP")
  )

  scores <- calc_sc_gsea_score(deTable)

  # Very small p-values should result in large scores (more significant)
  expect_true(abs(scores["VerySmallP"]) > abs(scores["LargeP"]))
  expect_true(abs(scores["TinyP"]) > abs(scores["MediumP"]))
})

test_that("calc_sc_gsea_score percentage weighting works correctly", {
  # Test that higher percentages result in higher absolute scores
  deTable <- data.frame(
    p_val = c(0.01, 0.01), # Same p-value
    avg_log2FC = c(2.0, 2.0), # Same fold change
    pct.1 = c(0.9, 0.3), # Different pct.1
    pct.2 = c(0.1, 0.2), # Different pct.2
    p_val_adj = c(0.05, 0.05),
    row.names = c("HighPct", "LowPct")
  )

  scores <- calc_sc_gsea_score(deTable)

  # HighPct should have higher score due to max(0.9, 0.1) = 0.9 vs max(0.3, 0.2) = 0.3
  expect_true(scores["HighPct"] > scores["LowPct"])
})

test_that("calc_sc_gsea_score sign handling works correctly", {
  # Test positive and negative fold changes
  deTable <- data.frame(
    p_val = c(0.01, 0.01),
    avg_log2FC = c(2.0, -2.0), # Opposite signs
    pct.1 = c(0.8, 0.8),
    pct.2 = c(0.2, 0.2),
    p_val_adj = c(0.05, 0.05),
    row.names = c("Positive", "Negative")
  )

  scores <- calc_sc_gsea_score(deTable)

  # Positive fold change should result in positive score component from sign
  expect_true(scores["Positive"] > 0)
  # Negative fold change should result in negative score component from sign
  expect_true(scores["Negative"] < 0)
  # Magnitudes should be similar (same p-val and pct values)
  expect_equal(
    abs(as.numeric(scores["Positive"])),
    abs(as.numeric(scores["Negative"])),
    tolerance = 1.0
  )
})
