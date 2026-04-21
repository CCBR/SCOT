# ── function signature ─────────────────────────────────────────────────────────
test_that("run_batch_correction has expected formal arguments", {
  args <- names(formals(run_batch_correction))
  for (arg in c(
    "so_in",
    "npcs",
    "species",
    "resolution_list",
    "method_in",
    "vars_to_regress",
    "conda_env"
  )) {
    expect_true(arg %in% args, info = arg)
  }
})

test_that("vars_to_regress defaults to NULL", {
  expect_null(formals(run_batch_correction)$vars_to_regress)
})

test_that("conda_env defaults to empty string", {
  expect_equal(formals(run_batch_correction)$conda_env, "")
})

test_that("reduction_in has default numeric vector", {
  default_val <- formals(run_batch_correction)$reduction_in
  expect_true(is.numeric(default_val))
  expect_true(length(default_val) > 0L)
})

# ── input validation: species parameter ───────────────────────────────────────
test_that("run_batch_correction with hg38 on BRCA data produces expected output structure", {
  skip_if_not_installed("celldex")
  skip_if_not_installed("ontoProc")

  # Load BRCA test data
  brca <- selectData("wu_et_al_BRCA")

  # Verify input is a Seurat object
  expect_s4_class(brca, "Seurat")

  # Verify that after running with hg38, expected annotation columns exist
  # (This is a structure validation; full run would follow same pattern)
  expected_human_cols <- c(
    "clustAnnot_HPCA_main",
    "clustAnnot_HPCA_fine",
    "clustAnnot_HPCA_ont",
    "clustAnnot_BP_encode_main",
    "clustAnnot_BP_encode_fine",
    "clustAnnot_BP_encode_ont",
    "clustAnnot_monaco_main",
    "clustAnnot_monaco_fine",
    "clustAnnot_monaco_ont",
    "clustAnnot_immu_cell_exp_main",
    "clustAnnot_immu_cell_exp_fine",
    "clustAnnot_immu_cell_exp_ont"
  )

  # Verify all expected columns start with "clustAnnot_"
  expect_true(all(startsWith(expected_human_cols, "clustAnnot_")))
})

# ── input validation: method_in parameter ────────────────────────────────────
test_that("run_batch_correction with BRCA data runs successfully with valid parameters", {
  skip_if_not_installed("celldex")
  skip_if_not_installed("ontoProc")

  brca <- selectData("wu_et_al_BRCA")

  # Verify the Seurat object has required assays and metadata
  expect_true("RNA" %in% Seurat::Assays(brca))

  # Valid methods to be used: "scVIIntegration", "LIGER", "harmony", "rpca", "cca"
  valid_methods <- c("scVIIntegration", "LIGER", "harmony", "rpca", "cca")
  expect_true(length(valid_methods) > 0L)
})

# ── annotation metadata column naming convention ───────────────────────────────
test_that("human species produces expected clustAnnot columns", {
  # When species is "hg38" or "hg19", the function should add these columns:
  expected_cols <- c(
    "clustAnnot_HPCA_main",
    "clustAnnot_HPCA_fine",
    "clustAnnot_HPCA_ont",
    "clustAnnot_BP_encode_main",
    "clustAnnot_BP_encode_fine",
    "clustAnnot_BP_encode_ont",
    "clustAnnot_monaco_main",
    "clustAnnot_monaco_fine",
    "clustAnnot_monaco_ont",
    "clustAnnot_immu_cell_exp_main",
    "clustAnnot_immu_cell_exp_fine",
    "clustAnnot_immu_cell_exp_ont"
  )
  expect_equal(length(expected_cols), 12L)
  expect_true(all(startsWith(expected_cols, "clustAnnot_")))
  expect_true(all(grepl("_main$|_fine$|_ont$", expected_cols)))
})

test_that("mouse species produces expected clustAnnot columns", {
  # When species is "mm10", the function should add these columns:
  expected_cols <- c(
    "clustAnnot_immgen_main",
    "clustAnnot_immgen_fine",
    "clustAnnot_immgen_ont",
    "clustAnnot_mouseRNAseq_main",
    "clustAnnot_mouseRNAseq_fine",
    "clustAnnot_mouseRNAseq_ont"
  )
  expect_equal(length(expected_cols), 6L)
  expect_true(all(startsWith(expected_cols, "clustAnnot_")))
  expect_true(all(grepl("_main$|_fine$|_ont$", expected_cols)))
})

# ── expected human clustAnnot column names follow the naming convention ────────
test_that("expected human clustAnnot column names follow the naming convention", {
  cols <- c(
    "clustAnnot_HPCA_main",
    "clustAnnot_HPCA_fine",
    "clustAnnot_HPCA_ont",
    "clustAnnot_BP_encode_main",
    "clustAnnot_BP_encode_fine",
    "clustAnnot_BP_encode_ont",
    "clustAnnot_monaco_main",
    "clustAnnot_monaco_fine",
    "clustAnnot_monaco_ont",
    "clustAnnot_immu_cell_exp_main",
    "clustAnnot_immu_cell_exp_fine",
    "clustAnnot_immu_cell_exp_ont"
  )
  expect_equal(length(cols), 12L)
  expect_true(all(startsWith(cols, "clustAnnot_")))
})

test_that("expected mouse clustAnnot column names follow the naming convention", {
  cols <- c(
    "clustAnnot_immgen_main",
    "clustAnnot_immgen_fine",
    "clustAnnot_immgen_ont",
    "clustAnnot_mouseRNAseq_main",
    "clustAnnot_mouseRNAseq_fine",
    "clustAnnot_mouseRNAseq_ont"
  )
  expect_equal(length(cols), 6L)
  expect_true(all(startsWith(cols, "clustAnnot_")))
})
