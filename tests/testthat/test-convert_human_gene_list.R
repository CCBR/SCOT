test_that("convert_human_gene_list returns character vector", {
  converted_result <- convert_human_gene_list(c("TP53", "EGFR", "CD3D"))
  expect_type(converted_result, "character")
})

test_that("known human genes map to expected mouse symbols", {
  converted_result <- convert_human_gene_list(c("TP53", "EGFR"))
  expect_true("Trp53" %in% converted_result)
  expect_true("Egfr" %in% converted_result)
})

test_that("non-character input raises an error", {
  expect_error(convert_human_gene_list(123), "character|genes")
  expect_error(convert_human_gene_list(list("TP53")), "character|genes")
  expect_error(convert_human_gene_list(NULL), "character|genes")
})

test_that("empty character vector returns empty character vector", {
  converted_result <- convert_human_gene_list(character(0))
  expect_identical(converted_result, character(0))
})

test_that("unrecognised gene symbols raise a valid-keys error", {
  expect_error(
    convert_human_gene_list(c("NOT_A_GENE")),
    "None of the keys entered are valid keys for 'SYMBOL'"
  )
})

test_that("mixed known and unknown genes return character vector with Trp53 present", {
  converted_result <- convert_human_gene_list(c("TP53", "NOT_A_GENE"))
  expect_type(converted_result, "character")
  expect_true("Trp53" %in% converted_result)
})

test_that("duplicate inputs produce same output", {
  genes <- c("TP53", "TP53", "EGFR")
  expect_identical(
    convert_human_gene_list(genes),
    convert_human_gene_list(genes)
  )
})

test_that("human genes from BRCA fixture convert to mouse title-case symbols", {
  brca_genes <- c(
    "TP53",
    "EGFR",
    "BRCA1",
    "BRCA2",
    "MED21",
    "TMPRSS5",
    "GPR139",
    "RHBDL3",
    "NDUFS8",
    "KLC1",
    "LRRC31",
    "GCFC2",
    "TCF19",
    "ASIC3",
    "STK32A",
    "PLCB1",
    "NOTCH1",
    "CDKN1A",
    "PTEN",
    "KRAS"
  )

  converted_result <- convert_human_gene_list(brca_genes)

  expect_type(converted_result, "character")

  # Genes starting with a letter must be in title case
  alpha_genes <- converted_result[grepl("^[A-Za-z]", converted_result)]
  expect_true(all(alpha_genes == stringr::str_to_title(alpha_genes)))

  # Genes starting with a digit must be all lowercase
  num_genes <- converted_result[grepl("^[0-9]", converted_result)]
  if (length(num_genes) > 0) {
    expect_true(all(grepl("^[0-9][a-z0-9]*$", num_genes)))
  }
})
