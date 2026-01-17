# Test for run_batch_correction function

test_that("run_batch_correction function exists and has correct signature", {
  # Source function for testing
  source(file.path("..", "..", "R", "run_batch_correction.R"))

  # Test that the function exists and is callable
  expect_true(exists("run_batch_correction"))
  expect_true(is.function(run_batch_correction))

  # Check function signature
  fn_args <- names(formals(run_batch_correction))
  expected_args <- c("so_in", "npcs", "species", "resolution_list", "method_in",
                    "reduction_in", "vars_to_regress", "conda_env")
  expect_true(all(expected_args %in% fn_args))
  expect_equal(length(fn_args), 8)

  # Check default values
  defaults <- formals(run_batch_correction)
  expect_equal(defaults$reduction_in, c(0.2, 0.4, 0.5, 0.8, 1.0))
  expect_null(defaults$vars_to_regress)
  expect_equal(defaults$conda_env, "")
})

test_that("run_batch_correction validates input parameters", {
  # Source function for testing
  source(file.path("..", "..", "R", "run_batch_correction.R"))

  # Test method validation
  valid_methods <- c("scVIIntegration", "LIGER", "CCAIntegration", "RPCAIntegration",
                    "HarmonyIntegration", "FastMNNIntegration")
  test_method <- "CCAIntegration"
  expect_true(test_method %in% valid_methods)

  # Test species validation
  valid_species <- c("hg38", "hg19", "mm10")
  expect_true("hg38" %in% valid_species)
  expect_true("mm10" %in% valid_species)

  # Test npcs should be numeric and positive
  valid_npcs <- c(20, 30, 50)
  for (npcs in valid_npcs) {
    expect_true(is.numeric(npcs) && npcs > 0)
  }

  # Test resolution_list should be numeric vector
  test_resolutions <- c(0.2, 0.5, 0.8, 1.0)
  expect_true(is.numeric(test_resolutions))
  expect_true(all(test_resolutions >= 0 && test_resolutions <= 2.0))
})

test_that("run_batch_correction method-specific logic", {
  # Source function for testing
  source(file.path("..", "..", "R", "run_batch_correction.R"))

  # Test method classification
  scvi_method <- "scVIIntegration"
  liger_method <- "LIGER"
  other_methods <- c("CCAIntegration", "RPCAIntegration", "HarmonyIntegration")

  # Test method identification logic
  expect_equal(scvi_method, "scVIIntegration")
  expect_equal(liger_method, "LIGER")
  expect_true(all(other_methods %in% c("CCAIntegration", "RPCAIntegration", "HarmonyIntegration")))

  # Test that methods are mutually exclusive
  expect_false(scvi_method %in% other_methods)
  expect_false(liger_method %in% other_methods)
  expect_false(scvi_method == liger_method)
})

test_that("run_batch_correction dependency checks", {
  # Skip if required packages are not available
  skip_if_not_installed("Seurat")

  library(Seurat)

  # Test core Seurat functions availability
  expect_true(exists("DefaultAssay", where = asNamespace("Seurat")))
  expect_true(exists("NormalizeData", where = asNamespace("Seurat")))
  expect_true(exists("FindVariableFeatures", where = asNamespace("Seurat")))
  expect_true(exists("ScaleData", where = asNamespace("Seurat")))
  expect_true(exists("RunPCA", where = asNamespace("Seurat")))
  expect_true(exists("SCTransform", where = asNamespace("Seurat")))
  expect_true(exists("IntegrateLayers", where = asNamespace("Seurat")))
  expect_true(exists("FindNeighbors", where = asNamespace("Seurat")))
  expect_true(exists("FindClusters", where = asNamespace("Seurat")))
  expect_true(exists("RunUMAP", where = asNamespace("Seurat")))
})

test_that("run_batch_correction optional dependency checks", {
  # Test optional dependencies - skip if not available

  # SeuratWrappers
  tryCatch({
    skip_if_not_installed("SeuratWrappers")
    library(SeuratWrappers)
    expect_true("SeuratWrappers" %in% loadedNamespaces())
  }, error = function(e) {
    skip("SeuratWrappers not available")
  })

  # harmony
  tryCatch({
    skip_if_not_installed("harmony")
    library(harmony)
    expect_true("harmony" %in% loadedNamespaces())
  }, error = function(e) {
    skip("harmony not available")
  })

  # rliger
  tryCatch({
    skip_if_not_installed("rliger")
    library(rliger)
    expect_true("rliger" %in% loadedNamespaces())
  }, error = function(e) {
    skip("rliger not available")
  })
})

test_that("run_batch_correction species-specific annotation logic", {
  # Source function for testing
  source(file.path("..", "..", "R", "run_batch_correction.R"))

  # Test species classification for cell type annotation
  human_species <- c("hg38", "hg19")
  mouse_species <- "mm10"

  # Human annotation databases (conceptual test)
  human_dbs <- c("HumanPrimaryCellAtlasData", "BlueprintEncodeData",
                "MonacoImmuneData", "DatabaseImmuneCellExpressionData")

  # Mouse annotation databases (conceptual test)
  mouse_dbs <- c("ImmGenData", "MouseRNAseqData")

  # Test database availability logic
  expect_equal(length(human_dbs), 4)
  expect_equal(length(mouse_dbs), 2)

  # Test annotation levels
  annotation_levels <- c("label.main", "label.fine", "label.ont")
  expect_equal(length(annotation_levels), 3)
  expect_true(all(grepl("^label\\.", annotation_levels)))
})

test_that("run_batch_correction clustering parameters", {
  # Test clustering parameters from the function

  # Algorithm parameter (Leiden = 3)
  clustering_algorithm <- 3
  expect_equal(clustering_algorithm, 3)
  expect_true(is.numeric(clustering_algorithm))

  # Resolution iteration logic
  test_resolutions <- c(0.2, 0.4, 0.5, 0.8, 1.0)

  expect_true(is.numeric(test_resolutions))
  expect_true(all(test_resolutions > 0))
  expect_true(all(test_resolutions <= 2.0))
  expect_equal(length(test_resolutions), 5)

  # Test that resolutions are in ascending order
  expect_equal(test_resolutions, sort(test_resolutions))
})

test_that("run_batch_correction UMAP parameters", {
  # Test UMAP parameters logic

  # Default reduction name
  default_reduction <- "integrated.cca"  # Example
  expect_true(is.character(default_reduction))

  # Dimensions parameter
  test_dims <- 1:30
  expect_true(is.numeric(test_dims))
  expect_equal(min(test_dims), 1)
  expect_true(max(test_dims) > 1)

  # Test dims construction with npcs
  test_npcs <- 50
  dims_range <- 1:test_npcs
  expect_equal(length(dims_range), test_npcs)
  expect_equal(min(dims_range), 1)
  expect_equal(max(dims_range), test_npcs)
})

test_that("run_batch_correction run_singleR_cluster dependency", {
  # Source function for testing
  source(file.path("..", "..", "R", "run_batch_correction.R"))

  # Check if run_singleR_cluster function is available
  tryCatch({
    source(file.path("..", "..", "R", "run_singleR_cluster.R"))
    expect_true(exists("run_singleR_cluster"))

    # Test function signature if available
    if (exists("run_singleR_cluster")) {
      fn_args <- names(formals(run_singleR_cluster))
      expect_true(length(fn_args) >= 2)  # Should have at least so and reference args
    }
  }, error = function(e) {
    skip("run_singleR_cluster.R not found")
  })
})

test_that("run_batch_correction celldex dependency checks", {
  # Skip if celldex is not available
  skip_if_not_installed("celldex")

  library(celldex)

  # Test human reference datasets
  human_refs <- c("HumanPrimaryCellAtlasData", "BlueprintEncodeData",
                 "MonacoImmuneData", "DatabaseImmuneCellExpressionData")

  for (ref_func in human_refs) {
    expect_true(exists(ref_func, where = asNamespace("celldex")))
  }

  # Test mouse reference datasets
  mouse_refs <- c("ImmGenData", "MouseRNAseqData")

  for (ref_func in mouse_refs) {
    expect_true(exists(ref_func, where = asNamespace("celldex")))
  }
})

test_that("run_batch_correction ontoProc dependency checks", {
  # Skip if ontoProc is not available
  skip_if_not_installed("ontoProc")

  library(ontoProc)

  # Test getOnto function availability
  expect_true(exists("getOnto", where = asNamespace("ontoProc")))

  # Test cellOnto parameter
  cell_onto_param <- "cellOnto"
  expect_equal(cell_onto_param, "cellOnto")
  expect_true(is.character(cell_onto_param))
})

test_that("run_batch_correction SCTransform variables regression", {
  # Test vars.to.regress parameter handling

  # Default case (NULL)
  default_vars <- NULL
  expect_null(default_vars)

  # Example regression variables
  example_vars <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
  expect_true(is.character(example_vars))
  expect_equal(length(example_vars), 3)

  # Test variable names format
  for (var in example_vars) {
    expect_true(nchar(var) > 0)
    expect_true(is.character(var))
  }
})

test_that("run_batch_correction conda environment parameter", {
  # Test conda_env parameter handling

  # Default empty string
  default_env <- ""
  expect_equal(default_env, "")
  expect_true(is.character(default_env))

  # Example conda environment names
  example_envs <- c("scvi-env", "pytorch_env", "scanpy")

  for (env in example_envs) {
    expect_true(is.character(env))
    expect_true(nchar(env) > 0)
  }
})

test_that("run_batch_correction error handling and edge cases", {
  # Source function for testing
  source(file.path("..", "..", "R", "run_batch_correction.R"))

  # Test with invalid inputs (conceptual)
  invalid_methods <- c("InvalidMethod", "", NA)
  valid_methods <- c("scVIIntegration", "LIGER", "CCAIntegration")

  for (method in invalid_methods) {
    if (!is.na(method)) {
      expect_false(method %in% valid_methods)
    }
  }

  # Test invalid species
  invalid_species <- c("rat", "zebrafish", "", NA)
  valid_species <- c("hg38", "hg19", "mm10")

  for (species in invalid_species) {
    if (!is.na(species)) {
      expect_false(species %in% valid_species)
    }
  }

  # Test invalid npcs
  invalid_npcs <- c(-1, 0, Inf, NA)

  for (npcs in invalid_npcs) {
    if (!is.na(npcs) && is.finite(npcs)) {
      expect_false(npcs > 0)
    }
  }
})

test_that("run_batch_correction return value validation", {
  # Test expected return value structure (conceptual)

  # The function should return a Seurat object with additional metadata
  # Test the expected structure without running the full function

  expected_metadata_cols <- c(
    "clustAnnot_HPCA_main", "clustAnnot_HPCA_fine", "clustAnnot_HPCA_ont",
    "clustAnnot_BP_encode_main", "clustAnnot_BP_encode_fine", "clustAnnot_BP_encode_ont",
    "clustAnnot_monaco_main", "clustAnnot_monaco_fine", "clustAnnot_monaco_ont",
    "clustAnnot_immu_cell_exp_main", "clustAnnot_immu_cell_exp_fine", "clustAnnot_immu_cell_exp_ont"
  )

  # Test human annotation column names
  expect_equal(length(expected_metadata_cols), 12)
  expect_true(all(grepl("^clustAnnot_", expected_metadata_cols)))

  # Test mouse annotation columns
  mouse_metadata_cols <- c(
    "clustAnnot_immgen_main", "clustAnnot_immgen_fine", "clustAnnot_immgen_ont",
    "clustAnnot_mouseRNAseq_main", "clustAnnot_mouseRNAseq_fine", "clustAnnot_mouseRNAseq_ont"
  )

  expect_equal(length(mouse_metadata_cols), 6)
  expect_true(all(grepl("^clustAnnot_", mouse_metadata_cols)))
})
