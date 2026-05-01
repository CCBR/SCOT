# Test for calc_sc_gsea_score function

test_that("calc_sc_gsea_score computes scores correctly", {
  # Create a mock differential expression table
  de_table <- data.frame(
    p_val = c(0.001, 0.05, 0.1, 0.0001),
    avg_log2FC = c(2.5, -1.2, 0.8, -3.1),
    pct.1 = c(0.8, 0.3, 0.6, 0.9),
    pct.2 = c(0.2, 0.7, 0.4, 0.1),
    p_val_adj = c(0.01, 0.15, 0.25, 0.001),
    row.names = c("Gene1", "Gene2", "Gene3", "Gene4")
  )

  # Calculate scores
  gsea_score_vect <- calc_sc_gsea_score(de_table)

  # Test basic properties
  expect_type(gsea_score_vect, "double")
  expect_equal(length(gsea_score_vect), nrow(de_table))
  expect_equal(names(gsea_score_vect), rownames(de_table))

  expect_equal(
    gsea_score_vect,
    c(
      Gene1 = 4.9,
      Gene2 = -2.11072099696479,
      Gene3 = 1.4,
      Gene4 = -6.7
    )
  )
})

test_that("calc_sc_gsea_score handles zero p-values correctly", {
  # Create table with zero p-values
  de_table <- data.frame(
    p_val = c(0, 0.05, 0),
    avg_log2FC = c(2.0, -1.5, -2.5),
    pct.1 = c(0.9, 0.4, 0.7),
    pct.2 = c(0.1, 0.6, 0.3),
    p_val_adj = c(0, 0.1, 0),
    row.names = c("ZeroGene1", "NonZeroGene", "ZeroGene2")
  )

  gsea_score_vect <- calc_sc_gsea_score(de_table)

  expect_equal(
    gsea_score_vect,
    c(
      ZeroGene1 = 452,
      NonZeroGene = -2.28061799739839,
      ZeroGene2 = -352.5
    )
  )
})

test_that("calc_sc_gsea_score handles edge cases", {
  # Test with minimal data - single gene
  de_table_minimal <- data.frame(
    p_val = 0.5,
    avg_log2FC = 1.0,
    pct.1 = 0.5,
    pct.2 = 0.3,
    p_val_adj = 0.6,
    row.names = "SingleGene"
  )

  gsea_score_vect_minimal <- calc_sc_gsea_score(de_table_minimal)
  expect_equal(gsea_score_vect_minimal, c(SingleGene = 1.15051499783199))

  # Test with zero fold change
  de_table_zero_fc <- data.frame(
    p_val = c(0.01, 0.05),
    avg_log2FC = c(0, 0),
    pct.1 = c(0.6, 0.4),
    pct.2 = c(0.4, 0.6),
    p_val_adj = c(0.05, 0.1),
    row.names = c("ZeroFC1", "ZeroFC2")
  )

  gsea_score_vect_zero_fc <- calc_sc_gsea_score(de_table_zero_fc)
  expect_equal(gsea_score_vect_zero_fc, c(ZeroFC1 = 0, ZeroFC2 = 0))

  # Test with NA values - should result in NA scores
  de_table_na <- data.frame(
    p_val = c(0.01, NA),
    avg_log2FC = c(1.5, 2.0),
    pct.1 = c(0.7, 0.8),
    pct.2 = c(0.3, 0.2),
    p_val_adj = c(0.05, NA),
    row.names = c("Gene1", "NA_Gene2")
  )

  gsea_score_vect_na <- calc_sc_gsea_score(de_table_na)
  expect_equal(gsea_score_vect_na, c(Gene1 = 2.9, NA_Gene2 = NA_real_))
})

test_that("calc_sc_gsea_score validates input structure", {
  # Test with missing required columns (missing pct.1, pct.2, p_val_adj)
  incomplete_table <- data.frame(
    p_val = c(0.01, 0.05),
    avg_log2FC = c(1.0, -0.5),
    row.names = c("Gene1", "Gene2")
  )

  expect_error(
    calc_sc_gsea_score(incomplete_table),
    regexp = "de_table is missing required columns"
  )

  # Test with NULL input (not a data.frame)
  expect_error(
    calc_sc_gsea_score(NULL),
    regexp = "de_table is not a data.frame"
  )

  # Test with non-data.frame input (string)
  expect_error(
    calc_sc_gsea_score("not_a_dataframe"),
    regexp = "de_table is not a data.frame"
  )

  # Test with non-data.frame input (vector)
  expect_error(
    calc_sc_gsea_score(c(1, 2, 3)),
    regexp = "de_table is not a data.frame"
  )

  # Test with missing single required column
  missing_pct2 <- data.frame(
    p_val = c(0.01, 0.05),
    avg_log2FC = c(1.0, -0.5),
    pct.1 = c(0.3, 0.5),
    p_val_adj = c(0.05, 0.10),
    row.names = c("Gene1", "Gene2")
  )
  expect_error(
    calc_sc_gsea_score(missing_pct2),
    regexp = "de_table is missing required columns.*pct.2"
  )
})

test_that("calc_sc_gsea_score handles various p-value ranges", {
  # Test with very small and large p-values
  de_table <- data.frame(
    p_val = c(1e-10, 0.99, 1e-300, 0.5),
    avg_log2FC = c(3.0, -1.0, 2.5, -0.5),
    pct.1 = c(0.8, 0.2, 0.9, 0.4),
    pct.2 = c(0.1, 0.8, 0.1, 0.6),
    p_val_adj = c(1e-8, 0.99, 1e-250, 0.7),
    row.names = c("VerySmallP", "LargeP", "TinyP", "MediumP")
  )

  gsea_score_vect <- calc_sc_gsea_score(de_table)

  expect_equal(
    gsea_score_vect,
    c(
      VerySmallP = 11,
      LargeP = -1.00349184432578,
      TinyP = 272.5,
      MediumP = -0.680617997398389
    )
  )
})

test_that("calc_sc_gsea_score percentage weighting works correctly", {
  # Test that higher percentages result in higher absolute scores
  de_table <- data.frame(
    p_val = c(0.01, 0.01), # Same p-value
    avg_log2FC = c(2.0, 2.0), # Same fold change
    pct.1 = c(0.9, 0.3), # Different pct.1
    pct.2 = c(0.1, 0.2), # Different pct.2
    p_val_adj = c(0.05, 0.05),
    row.names = c("HighPct", "LowPct")
  )

  gsea_score_vect <- calc_sc_gsea_score(de_table)

  expect_equal(gsea_score_vect, c(HighPct = 3.8, LowPct = 2.6))
})

test_that("calc_sc_gsea_score sign handling works correctly", {
  # Test positive and negative fold changes
  de_table <- data.frame(
    p_val = c(0.01, 0.01),
    avg_log2FC = c(2.0, -2.0), # Opposite signs
    pct.1 = c(0.8, 0.8),
    pct.2 = c(0.2, 0.2),
    p_val_adj = c(0.05, 0.05),
    row.names = c("Positive", "Negative")
  )

  gsea_score_vect <- calc_sc_gsea_score(de_table)

  expect_equal(gsea_score_vect, c(Positive = 3.6, Negative = -3.6))
})

test_that("calc_sc_gsea_score preserves provided row names", {
  de_table_empty_rowname <- data.frame(
    p_val = c(0.01, 0.05),
    avg_log2FC = c(1.0, -1.0),
    pct.1 = c(0.5, 0.4),
    pct.2 = c(0.2, 0.6),
    p_val_adj = c(0.05, 0.1),
    row.names = c("", "Gene2")
  )

  gsea_score_vect_empty <- calc_sc_gsea_score(de_table_empty_rowname)
  expect_equal(names(gsea_score_vect_empty), c("", "Gene2"))

  de_table_dup_rowname <- data.frame(
    p_val = c(0.02, 0.03),
    avg_log2FC = c(0.9, -1.2),
    pct.1 = c(0.7, 0.4),
    pct.2 = c(0.3, 0.6),
    p_val_adj = c(0.08, 0.1)
  )
  attr(de_table_dup_rowname, "row.names") <- c("DupGene", "DupGene")

  gsea_score_vect_dup <- calc_sc_gsea_score(de_table_dup_rowname)
  expect_equal(names(gsea_score_vect_dup), c("DupGene", "DupGene"))
})
