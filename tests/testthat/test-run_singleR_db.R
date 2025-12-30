# Tests for run_singleR_db

context("run_singleR_db")

test_that("run_singleR_db exists and has correct signature", {
  source(file.path("..", "..", "R", "run_singleR_db.R"))
  expect_true(exists("run_singleR_db"))
  expect_true(is.function(run_singleR_db))
  args <- names(formals(run_singleR_db))
  expect_equal(args, c("so_in", "species"))
})

test_that("run_singleR_db dependencies available", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("celldex")
  skip_if_not_installed("ontoProc")
  skip_if_not_installed("SingleR")
})

make_mock_seurat_min <- function(n_genes = 20, n_cells = 12) {
  library(Seurat)
  library(Matrix)
  set.seed(123)
  L <- n_genes * 2
  counts <- Matrix::sparseMatrix(
    i = sample(seq_len(n_genes), L, replace = TRUE),
    j = sample(seq_len(n_cells), L, replace = TRUE),
    x = rpois(L, 2),
    dims = c(n_genes, n_cells)
  )
  rownames(counts) <- paste0("Gene", seq_len(n_genes))
  colnames(counts) <- paste0("Cell", seq_len(n_cells))
  so <- Seurat::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
  Seurat::DefaultAssay(so) <- "RNA"
  so <- Seurat::NormalizeData(so)
  so <- Seurat::FindVariableFeatures(so)
  so <- Seurat::ScaleData(so)
  so <- Seurat::RunPCA(so, npcs = 10)
  so
}

test_that("run_singleR_db species parameter validation", {
  source(file.path("..", "..", "R", "run_singleR_db.R"))
  expect_true(all(c("hg38", "hg19", "mm10") %in% c("hg38", "hg19", "mm10")))
  expect_false("human" %in% c("hg38", "hg19", "mm10"))
})

test_that("run_singleR_db conceptual annotation columns", {
  source(file.path("..", "..", "R", "run_singleR_db.R"))
  # Human columns
  human_cols <- c(
    "HPCA_main", "HPCA_fine", "HPCA_ont",
    "BP_encode_main", "BP_encode_fine", "BP_encode_ont",
    "monaco_main", "monaco_fine", "monaco_ont",
    "immu_cell_exp_main", "immu_cell_exp_fine", "immu_cell_exp_ont",
    "annot"
  )
  expect_true(all(grepl("^HPCA_|^BP_encode_|^monaco_|^immu_cell_exp_|^annot$", human_cols)))
  
  # Mouse columns
  mouse_cols <- c(
    "immgen_main", "immgen_fine", "immgen_ont",
    "mouseRNAseq_main", "mouseRNAseq_fine", "mouseRNAseq_ont",
    "annot"
  )
  expect_true(all(grepl("^immgen_|^mouseRNAseq_|^annot$", mouse_cols)))
})

test_that("run_singleR_db errors gracefully without SCT/refs (conceptual)", {
  skip_if_not_installed("Seurat")
  source(file.path("..", "..", "R", "run_singleR_db.R"))
  so <- make_mock_seurat_min(20, 12)
  expect_error(run_singleR_db(so, "hg38"))
})
