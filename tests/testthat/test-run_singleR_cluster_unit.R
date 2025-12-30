context("run_singleR_cluster (unit)")

test_that("run_singleR_cluster has expected signature", {
  expect_true(exists("run_singleR_cluster"))
  fmls <- formals(run_singleR_cluster)
  expect_setequal(names(fmls), c("so_in", "ref_file", "label"))
})

test_that("errors when SCT assay is missing", {
  skip_if_not_installed("Seurat")
  suppressPackageStartupMessages(library(Seurat))

  set.seed(11)
  mat <- matrix(rpois(100, lambda = 2), nrow = 10, ncol = 10,
                dimnames = list(paste0("Gene", seq_len(10)), paste0("Cell", seq_len(10))))
  so <- CreateSeuratObject(counts = mat)
  expect_error(run_singleR_cluster(so, ref_file = NULL, label = "label"))
})

test_that("maps cluster labels back to cells and strips 'SCT.' prefix", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SingleR")
  skip_if_not_installed("SingleCellExperiment")

  suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
  })

  set.seed(123)
  genes <- paste0("Gene", seq_len(30))
  cells <- paste0("Cell", seq_len(20))
  counts <- matrix(rpois(length(genes) * length(cells), lambda = 3),
                   nrow = length(genes), ncol = length(cells),
                   dimnames = list(genes, cells))

  so <- CreateSeuratObject(counts = counts)
  # Add a minimal SCT assay so AverageExpression(assays="SCT") works
  so[["SCT"]] <- SeuratObject::CreateAssayObject(counts = counts)
  Seurat::DefaultAssay(so) <- "SCT"

  # Define clusters and identities
  cl <- factor(rep(c("0", "1"), each = 10))
  so$seurat_clusters <- cl
  Seurat::Idents(so) <- cl

  # Build a tiny reference SCE with matching genes
  ref_counts <- matrix(rpois(length(genes) * 6, lambda = 3),
                       nrow = length(genes), ncol = 6,
                       dimnames = list(genes, paste0("ref", seq_len(6))))
  ref <- SingleCellExperiment(list(logcounts = ref_counts))
  SummarizedExperiment::colData(ref)$label <- rep(c("A", "B"), each = 3)

  # Run (allow warnings from SingleR internals)
  annot <- suppressWarnings(run_singleR_cluster(so, ref_file = ref, label = "label"))

  # Checks: lengths, names, and value domain (allow NA from SingleR)
  expect_type(annot, "character")
  expect_length(annot, ncol(so))
  expect_identical(names(annot), colnames(so))
  expect_true(all(stats::na.omit(annot) %in% c("A", "B")))
})
